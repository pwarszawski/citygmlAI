// This software is based on pugixml library (http://pugixml.org).
// pugixml is Copyright (C) 2006-2018 Arseny Kapoulkine.

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <map>
#include <set>
#include "pugixml.hpp"

namespace fs = std::filesystem;

// ----------------------------------------------------------------
// Konfiguracja pobierana z pliku .ini
// ----------------------------------------------------------------
struct Config {
    std::string trackFilename = "EXPORT.SCN";
    std::string inputDir = "CityGML-walbrzych";
    std::string outputFilename = "citygml.scm";
    double filterDistance = 2000.0; 
    bool swapGmlXy = true;        
    double minVertexDistSq = 0.01 * 0.01; 
    double minTriangleArea = 0.01;
    std::string texRoof = "roof/karpiowka";
    std::string texWall = "beton2";
    std::string texGround = "asphaltdark1";

    bool load(const std::string& iniPath) {
        std::ifstream file(iniPath);
        if (!file.is_open()) {
            std::cerr << "BLAD: Nie znaleziono pliku konfiguracyjnego '" << iniPath << "'!\n";
            std::cerr << "Zbuduj plik .ini z parametrami przed uruchomieniem programu.\n";
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            // Ignorowanie komentarzy (# lub ;)
            auto commentPos = line.find('#');
            if (commentPos == std::string::npos) commentPos = line.find(';');
            if (commentPos != std::string::npos) line = line.substr(0, commentPos);

            // Usuwanie bialych znakow
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            // Ignorowanie pustych linii i naglowkow sekcji
            if (line.empty() || line[0] == '[') continue;

            auto delim = line.find('=');
            if (delim != std::string::npos) {
                std::string key = line.substr(0, delim);
                std::string val = line.substr(delim + 1);
                
                key.erase(key.find_last_not_of(" \t") + 1);
                val.erase(0, val.find_first_not_of(" \t"));

                try {
                    if (key == "TRACK_FILENAME") trackFilename = val;
                    else if (key == "INPUT_DIR") inputDir = val;
                    else if (key == "OUTPUT_FILENAME") outputFilename = val;
                    else if (key == "FILTER_DISTANCE") filterDistance = std::stod(val);
                    else if (key == "SWAP_GML_XY") swapGmlXy = (val == "true" || val == "1");
                    else if (key == "MIN_VERTEX_DIST") { double v = std::stod(val); minVertexDistSq = v * v; }
                    else if (key == "MIN_TRIANGLE_AREA") minTriangleArea = std::stod(val);
                    else if (key == "TEX_ROOF") texRoof = val;
                    else if (key == "TEX_WALL") texWall = val;
                    else if (key == "TEX_GROUND") texGround = val;
                } catch(...) {
                    std::cerr << "Blad parsowania wartosci dla klucza: " << key << "\n";
                }
            }
        }
        return true;
    }

    void print() const {
        std::cout << "=========================================\n";
        std::cout << " citygmlAI v28\n";
        std::cout << "=========================================\n";
        std::cout << " Wczytana konfiguracja (citygmlAI.ini):\n";
        std::cout << "=========================================\n";
        std::cout << " -> Plik bazowy SCN : " << trackFilename << "\n";
        std::cout << " -> Katalog zrodlowy: " << inputDir << "\n";
        std::cout << " -> Plik wynikowy   : " << outputFilename << "\n";
        std::cout << " -> Max Dystans     : " << filterDistance << " m\n";
        std::cout << " -> Zamiana osi X/Y : " << (swapGmlXy ? "TAK" : "NIE") << "\n";
        std::cout << " -> Min Pow. Trojkata: " << minTriangleArea << " m2\n";
        std::cout << " -> Tekstura dachu  : " << texRoof << "\n";
        std::cout << " -> Tekstura sciany : " << texWall << "\n";
        std::cout << " -> Tekstura gruntu : " << texGround << "\n";
        std::cout << "=========================================\n\n";
    }
} cfg;

struct Point { double x, y, z; };
struct Vector3 { double nx, ny, nz; };

// ----------------------------------------------------------------
// Struktura przechowująca dane o powierzchniach do dalszej obróbki
// ----------------------------------------------------------------
struct PolyData {
    std::vector<Point> verts;
    Vector3 normal;
    std::string explicitType; // Typ tekstury wprost z XML (np. RoofSurface)
    double centerY;           // Środek wysokości do oceny czy to dach czy podłoga
};

// ----------------------------------------------------------------
// Funkcje pomocnicze
// ----------------------------------------------------------------
double parseDouble(const std::string& str) {
    std::string s = str;
    std::replace(s.begin(), s.end(), ',', '.');
    try { return std::stod(s); } catch (...) { return 0.0; }
}

double distanceSquared2D(double x1, double z1, double x2, double z2) {
    return (x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2);
}

double distSq3D(Point p1, Point p2) {
    return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z);
}

// ----------------------------------------------------------------
// Algorytm Newella do obliczania normalnych dla wielokątów
// ----------------------------------------------------------------
Vector3 calculateNewellNormal(const std::vector<Point>& poly) {
    Vector3 n = {0, 0, 0};
    for (size_t i = 0; i < poly.size(); ++i) {
        size_t next = (i + 1) % poly.size();
        n.nx += (poly[i].y - poly[next].y) * (poly[i].z + poly[next].z);
        n.ny += (poly[i].z - poly[next].z) * (poly[i].x + poly[next].x);
        n.nz += (poly[i].x - poly[next].x) * (poly[i].y + poly[next].y);
    }
    double len = std::sqrt(n.nx*n.nx + n.ny*n.ny + n.nz*n.nz);
    if (len > 0.00001) { n.nx /= len; n.ny /= len; n.nz /= len; }
    else { n = {0, 1, 0}; } 
    return n;
}

// ----------------------------------------------------------------
// Triangulacja metodą "ear clipping"
// ----------------------------------------------------------------
double getSignedArea(const std::vector<Point>& poly) {
    double area = 0.0;
    for (size_t i = 0; i < poly.size(); i++) {
        size_t j = (i + 1) % poly.size();
        area += (poly[j].x - poly[i].x) * (poly[j].y + poly[i].y);
    }
    return area / 2.0;
}

bool isPointInTriangle(Point p, Point a, Point b, Point c) {
    auto sign = [](Point p1, Point p2, Point p3) {
        return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
    };
    double d1 = sign(p, a, b);
    double d2 = sign(p, b, c);
    double d3 = sign(p, c, a);
    bool has_neg = (d1 < -0.001) || (d2 < -0.001) || (d3 < -0.001);
    bool has_pos = (d1 > 0.001) || (d2 > 0.001) || (d3 > 0.001);
    return !(has_neg && has_pos);
}

std::vector<int> triangulatePolygon(const std::vector<Point>& poly3d, Vector3 normal) {
    std::vector<int> indices;
    int n = poly3d.size();
    if (n < 3) return indices;

    std::vector<Point> poly2d(n);
    int dropAxis = 0; 
    if (std::abs(normal.ny) > std::abs(normal.nx) && std::abs(normal.ny) > std::abs(normal.nz)) dropAxis = 1;
    else if (std::abs(normal.nz) > std::abs(normal.nx)) dropAxis = 2;

    for(int i=0; i<n; ++i) {
        if (dropAxis == 0) { poly2d[i] = {poly3d[i].y, poly3d[i].z, 0}; }
        else if (dropAxis == 1) { poly2d[i] = {poly3d[i].x, poly3d[i].z, 0}; }
        else { poly2d[i] = {poly3d[i].x, poly3d[i].y, 0}; }
    }

    if (getSignedArea(poly2d) > 0) { 
        for(auto& p : poly2d) p.x = -p.x; 
    }

    std::vector<int> availPoints(n);
    for(int i=0; i<n; ++i) availPoints[i] = i;

    int count = n;
    int prev_count = 0;

    while (count > 2) {
        if (prev_count == count) break; 
        prev_count = count;

        for (int i = 0; i < count; ++i) {
            int i0 = availPoints[(i + count - 1) % count];
            int i1 = availPoints[i];
            int i2 = availPoints[(i + 1) % count];

            Point a = poly2d[i0];
            Point b = poly2d[i1];
            Point c = poly2d[i2];

            double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if (cross <= 0.00001) continue; 

            bool isEar = true;
            for (int j = 0; j < count; ++j) {
                if (j == (i + count - 1) % count || j == i || j == (i + 1) % count) continue;
                if (isPointInTriangle(poly2d[availPoints[j]], a, b, c)) {
                    isEar = false; break;
                }
            }

            if (isEar) {
                indices.push_back(i0); indices.push_back(i1); indices.push_back(i2);
                availPoints.erase(availPoints.begin() + i);
                count--;
                break;
            }
        }
    }
    
    // Jeśli triangulacja się nie powiodła (np. z powodu złożonego kształtu), zwróć prosty fan
    if (indices.empty() || count > 2) {
        indices.clear();
        for (int i = 1; i < n - 1; ++i) {
            indices.push_back(0); indices.push_back(i); indices.push_back(i + 1);
        }
    }
    return indices;
}

std::string getExplicitTexture(pugi::xml_node node) {
    pugi::xml_node current = node;
    for (int i = 0; i < 10; ++i) { 
        if (!current) break;
        std::string name = current.name();
        if (name.find("RoofSurface") != std::string::npos) return cfg.texRoof;
        if (name.find("GroundSurface") != std::string::npos) return cfg.texGround;
        if (name.find("WallSurface") != std::string::npos) return cfg.texWall;
        current = current.parent();
    }
    return ""; 
}

class CityGMLConverter {
    double offsetNorth = 0.0, offsetEast = 0.0;
    std::vector<Point> trackPoints; 
    struct BoundingBox {
        double minX = 1e18, maxX = -1e18, minZ = 1e18, maxZ = -1e18;
        void update(double x, double z) {
            if (x < minX) minX = x; if (x > maxX) maxX = x;
            if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
        }
        bool isInside(double x, double z, double r) const {
            return (x >= minX - r && x <= maxX + r && z >= minZ - r && z <= maxZ + r);
        }
    } trackBounds;

public:
    bool loadContextFromSCN(const std::string& scnPath) {
        std::ifstream file(scnPath);
        if (!file.is_open()) return false;
        std::string line; bool offsetFound = false, insideTrack = false;
        while (std::getline(file, line)) {
            if (!offsetFound && line.find("//$g") != std::string::npos) {
                std::stringstream ss(line); std::string t; std::vector<double> n;
                while (ss >> t) { try { n.push_back(std::stod(t)); } catch(...) {} }
                if (n.size() >= 2) { offsetEast = n[0] * 1000.0; offsetNorth = n[1] * 1000.0; offsetFound = true; }
                continue;
            }
            if (!insideTrack && line.find("node") == 0 && line.find("track") != std::string::npos) insideTrack = true;
            if (insideTrack && line.find("endtrack") == 0) insideTrack = false;
            if (insideTrack) {
                size_t f = line.find_first_not_of(" \t");
                if (f == std::string::npos || isalpha(line[f])) continue;
                std::stringstream ss(line); double x, y, z;
                if (ss >> x >> y >> z) { if (std::abs(x) > 0.1) { trackPoints.push_back({x, y, z}); trackBounds.update(x, z); } }
            }
        }
        std::cout << "   OFFSET WYZNACZONY: E=" << offsetEast << " N=" << offsetNorth << "\n";
        return offsetFound;
    }

    Point transform(double geoN, double geoE, double geoH) {
        return { offsetEast - geoE, geoH, geoN - offsetNorth };
    }

    void processFile(const std::string& gmlPath, std::ofstream& outFile, int& globalCounter) {
        pugi::xml_document doc;
        if (!doc.load_file(gmlPath.c_str())) return;
        std::cout << "Plik: " << fs::path(gmlPath).filename().string();

        auto buildings = doc.select_nodes("//*[local-name()='Building']");
        int added = 0;
        
        for (auto& bNode : buildings) {
            // ŚCIEŻKA A: LoD2 - PRIORYTET
            auto thematicNodes = bNode.node().select_nodes(".//*[local-name()='RoofSurface' or local-name()='WallSurface' or local-name()='GroundSurface' or local-name()='FloorSurface']");
            
            bool isLoD2 = (thematicNodes.size() > 0);
            std::vector<pugi::xml_node> surfacesToProcess;
            
            if (isLoD2) {
                for (auto& t : thematicNodes) surfacesToProcess.push_back(t.node());
            } else {
                // ŚCIEŻKA B: LoD1 - FALLBACK
                auto solidNodes = bNode.node().select_nodes(".//*[local-name()='posList']");
                for (auto& s : solidNodes) surfacesToProcess.push_back(s.node().parent());
            }

            if (surfacesToProcess.empty()) continue;

            std::vector<PolyData> polys;
            double buildMinY = 1e18; 
            Point firstPoint = {0,0,0}; 
            bool firstPointSet = false;

            // Etap 1: Parsowanie wierzchołków i obliczanie normalnych
            for (auto& surface : surfacesToProcess) {
                pugi::xml_node posNode = surface.child("gml:lod2MultiSurface");
                if (!posNode) posNode = surface;
                posNode = posNode.select_node(".//*[local-name()='posList']").node();
                if (!posNode) continue;

                std::stringstream ss(posNode.child_value());
                std::vector<double> raw; std::string t;
                while (ss >> t) raw.push_back(parseDouble(t));
                if (raw.size() < 9) continue;

                std::vector<Point> p_vec;
                double localMinY = 1e18;
                for(size_t i=0; i<raw.size(); i+=3) {
                    Point p;
                    if (cfg.swapGmlXy) p = transform(raw[i+1], raw[i], raw[i+2]);
                    else p = transform(raw[i], raw[i+1], raw[i+2]);
                    p_vec.push_back(p);
                    if (!firstPointSet) { firstPoint = p; firstPointSet = true; }
                    if (p.y < buildMinY) buildMinY = p.y; 
                    if (p.y < localMinY) localMinY = p.y;
                }

                if (p_vec.size() > 3 && 
                    std::abs(p_vec.front().x - p_vec.back().x) < 0.01 &&
                    std::abs(p_vec.front().z - p_vec.back().z) < 0.01) {
                    p_vec.pop_back();
                }
                if (p_vec.size() < 3) continue;

                Vector3 n = calculateNewellNormal(p_vec);

                PolyData pd;
                pd.verts = p_vec;
                pd.normal = n;
                pd.explicitType = getExplicitTexture(surface);
                pd.centerY = localMinY; 
                polys.push_back(pd);
            }

            if (polys.empty()) continue;

            if (!trackBounds.isInside(firstPoint.x, firstPoint.z, cfg.filterDistance)) continue;
            bool near = false;
            for(const auto& tp : trackPoints) {
                if (distanceSquared2D(firstPoint.x, firstPoint.z, tp.x, tp.z) < (cfg.filterDistance * cfg.filterDistance)) { 
                    near = true; break; 
                }
            }
            if (!near) continue;

            struct OutputMesh { std::vector<Point> v; std::vector<Vector3> n; };
            std::map<std::string, OutputMesh> groups;

            // Etap 2: Teksturowanie i triangulacja
            for (auto& pd : polys) {
                std::string tex = cfg.texWall;
                bool flipWinding = false;
                
                if (!pd.explicitType.empty()) {
                    tex = pd.explicitType;
                } else {
                    if (std::abs(pd.normal.ny) > 0.7) { 
                        if (pd.centerY <= buildMinY + 0.5) {
                            tex = cfg.texGround;
                        } else {
                            tex = cfg.texRoof;
                            if (pd.normal.ny < 0) flipWinding = true;
                        }
                    } else {
                        tex = cfg.texWall; 
                    }
                }
                
                if (tex == cfg.texRoof && pd.normal.ny < -0.1) flipWinding = true;
                
                if (flipWinding) {
                    pd.normal.nx = -pd.normal.nx; 
                    pd.normal.ny = -pd.normal.ny; 
                    pd.normal.nz = -pd.normal.nz;
                }

                std::vector<int> indices = triangulatePolygon(pd.verts, pd.normal);
                if (flipWinding) std::reverse(indices.begin(), indices.end());

                for (size_t i = 0; i < indices.size(); ++i) {
                    groups[tex].v.push_back(pd.verts[indices[i]]);
                    groups[tex].n.push_back(pd.normal);
                }
            }

            // Etap 3: Zapisywanie do pliku'
            for (auto& g : groups) {
                if (g.second.v.empty()) continue;
                
                outFile << "node -1 0 GML_" << globalCounter++ << " triangles " << g.first << "\n";
                size_t count = g.second.v.size();
                
                for (size_t i = 0; i < count; ++i) {
                    outFile << std::fixed << std::setprecision(2) 
                            << g.second.v[i].x << " " << g.second.v[i].y << " " << g.second.v[i].z << " " 
                            << std::setprecision(3) << g.second.n[i].nx << " " << g.second.n[i].ny << " " << g.second.n[i].nz << " 0 0";
                    
                    if (i < count - 1) {
                        outFile << " end";
                    }
                    outFile << "\n";
                }
                outFile << "endtri\n";
                added++;
            }
        }
        std::cout << " -> OK: " << added << "\n";
    }
};

// ----------------------------------------------------------------
// Main
// ----------------------------------------------------------------

int main() {
    if (!cfg.load("citygmlAI.ini")) {
        // Zatrzymanie pracy programu jeśli plik .ini nie zostal zaladowany
        return 1;
    }
    cfg.print();

    CityGMLConverter conv;
    if (!conv.loadContextFromSCN(cfg.trackFilename)) {
        std::cerr << "Nie udalo sie zaladowac pliku " << cfg.trackFilename << "\n";
        return 1;
    }

    std::ofstream outFile(cfg.outputFilename);
    outFile << "// Obiekty wygenerowano programem citygmlAI v28\n";
    
    int counter = 1;
    if (!fs::exists(cfg.inputDir)) {
        std::cerr << "Katalog " << cfg.inputDir << " nie istnieje!\n";
        return 1;
    }

    for (const auto& entry : fs::directory_iterator(cfg.inputDir)) {
        std::string ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".gml" || ext == ".xml") {
            conv.processFile(entry.path().string(), outFile, counter);
        }
    }
    
    std::cout << "\nZakonczono pomyslnie. Wygenerowano: " << counter - 1 << " obiektow.\n";
    return 0;
}