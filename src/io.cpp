#include "io.h"

//-- Input functions
bool IO::read_config(const char* file){ //todo
    return true;
}

bool IO::read_point_cloud(const char* file, Point_set_3& pc) {
    std::ifstream ifile(file, std::ios_base::binary);
    ifile >> pc;
    std::cerr << "POINT CLOUD: "<< pc.size() << " point read" << std::endl;
    return true;
}

bool IO::read_polygons(const char* file, nlohmann::json& j) {
    std::ifstream ifs(file);
    j = nlohmann::json::parse(ifs);
    return true;
}

//todo make output vectors depend on outputLayerID's rather than topoClasses - needed for many surfaceLayers
//-- Output functions
void IO::output_obj(const std::vector<TopoFeature*>& allFeatures) {
    using namespace config;
    int numOutputLayers = TopoFeature::get_num_output_layers();

    std::vector<std::ofstream> of;
    std::vector<std::string>   fs(numOutputLayers), bs(numOutputLayers);

    std::vector<std::unordered_map<std::string, unsigned long>> dPts(numOutputLayers);
    //-- Output points
    for (auto& f : allFeatures) {
        if (!f->is_active()) continue;
        if (outputSeparately)
            IO::get_obj_pts(f->get_mesh(),
                            fs[f->get_output_layer_id()],
                            bs[f->get_output_layer_id()],
                            dPts[f->get_output_layer_id()]);
        else
            IO::get_obj_pts(f->get_mesh(),
                            fs[f->get_output_layer_id()],
                            bs[f->get_output_layer_id()],
                            dPts[0]);
    }

    //-- Get output layer names -- TEMP, could add it to config and have user defined names for surf layers
    std::vector<std::string> outputLayerName = {"Terrain", "Buildings", "Sides", "Top"};
    int surfLayer = 1;
    for (auto i = 4; i < numOutputLayers; ++i) {
        outputLayerName.push_back("SurfaceLayer_" + std::to_string(surfLayer++));
    }

    //-- Add class name and output to file
    if (!outputSeparately) of.emplace_back().open(outputFileName + ".obj");
    for (int i = 0; i < fs.size(); ++i) {
        if (bs[i].empty()) continue;
        if (outputSeparately) of.emplace_back().open(outputFileName + "_" + outputLayerName[i] + ".obj");

        of.back() << fs[i] << "\ng " << outputLayerName[i] << bs[i];
    }
    for (auto& f : of) f.close();
}

void IO::output_stl(const std::vector<TopoFeature*>& allFeatures) {
    using namespace config;
    std::vector<std::ofstream> of;
    std::vector<std::string>   fs(topoClassName.size());

    //-- Get all triangles
    for (auto& f : allFeatures) {
        if (!f->is_active()) continue;
        IO::get_stl_pts(f->get_mesh(), fs[f->get_class()]);
    }

    //-- Add class name and output to file
    if (!outputSeparately) of.emplace_back().open(outputFileName + ".stl");
    for (int i = 0; i < fs.size(); ++i) {
        if (fs[i].empty()) continue;
        if (outputSeparately) of.emplace_back().open(outputFileName + "_" + topoClassName.at(i) + ".stl");

        of.back() << "\nsolid " << topoClassName.at(i);
        of.back() << fs[i];
        of.back() << "\nendsolid " << topoClassName.at(i);
    }
    for (auto& f : of) f.close();
}

void IO::output_cityjson(const std::vector<TopoFeature*>& allFeatures) {
    using namespace config;
    std::ofstream of;
    nlohmann::json j;

    j["type"] = "CityJSON";
    j["version"] = "1.0";
    j["metadata"] = {};
    std::vector<double> bbox = Boundary::get_domain_bbox();
    j["metadata"]["geographicalExtent"] = Boundary::get_domain_bbox();
    j["metadata"]["referenceSystem"] = "urn:ogc:def:crs:EPSG::7415";
    std::unordered_map<std::string, unsigned long> dPts;
    for (auto& f : allFeatures) {
        // Only Buildings and Terrain for now
        if (f->get_class() != BUILDING && f->get_class() != TERRAIN) continue;
        //-- Get feature info
        nlohmann::json b;
        f->get_cityjson_info(b);

        //-- Get feature geometry
        nlohmann::json g;
        IO::get_cityjson_geom(f->get_mesh(), g, dPts, f->get_cityjson_primitive());

        //-- Get feature semantics
        f->get_cityjson_semantics(g);

        //-- Append to main json struct
        b["geometry"].push_back(g);
        j["CityObjects"][f->get_id()] = b;
    }

    //-- Vertices - store them in a vector to quickly sort
    std::vector<std::string> thepts;
    thepts.resize(dPts.size());
    for (auto& p : dPts)
        thepts[p.second] = p.first;
    dPts.clear();
    for (auto& p : thepts) {
        std::vector<std::string> c;
        boost::split(c, p, boost::is_any_of(" "));
        j["vertices"].push_back({std::stod(c[0], NULL), std::stod(c[1], NULL), std::stod(c[2], NULL) });
    }

    of.open(outputFileName + ".json");
    of << j.dump() << std::endl;
}

void IO::get_obj_pts(const Mesh& mesh,
                     std::string& fs,
                     std::string& bs,
                     std::unordered_map<std::string, unsigned long>& dPts)
{
    for (auto& face : mesh.faces()) {
        std::vector<unsigned long> faceIdx; faceIdx.reserve(3);
        std::string fsTemp;
        std::string bsTemp;
        for (auto index : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            std::string pt = gen_key_bucket(mesh.point(index));
            auto it = dPts.find(pt);
            if (it == dPts.end()) {
                fs += "\nv " + pt;
                bsTemp += " " + std::to_string(dPts.size() + 1);

                faceIdx.push_back(dPts.size() + 1);

                dPts[pt] = dPts.size() + 1;
            } else {
                bsTemp += " " + std::to_string(it->second);

                faceIdx.push_back(it->second);
            }
        }
        //- Check for problematic faces
        std::sort(faceIdx.begin(), faceIdx.end());
        auto it = std::unique(faceIdx.begin(), faceIdx.end());
        bool wasUnique = (it == faceIdx.end());

        if (wasUnique) {
            bs += "\nf";
            bs += bsTemp;
        } else {
//            std::cerr << "Found duplicates!" << std::endl;
        }
    }
}

void IO::get_stl_pts(Mesh& mesh, std::string& fs) {
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    auto fnormals = mesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_normals(mesh, vnormals, fnormals);
    for (auto& face : mesh.faces()) {
        fs += "\nfacet normal " + gen_key_bucket(fnormals[face]);
        fs += "\n    outer loop";
        for (auto index : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            fs += "\n        vertex " + gen_key_bucket(mesh.point(index));
        }
        fs += "\n    endloop";
        fs += "\nendfacet";
    }
}

void IO::get_cityjson_geom(const Mesh& mesh, nlohmann::json& g, std::unordered_map<std::string, unsigned long>& dPts,
                           std::string primitive) {
    g["type"] = primitive;
    g["lod"] = config::lod;
    g["boundaries"];
    for (auto& face : mesh.faces()) {
        std::vector<unsigned long> faceIdx;  faceIdx.reserve(3);
        std::vector<unsigned long> tempPoly; tempPoly.reserve(3);
        for (auto index : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            std::string pt = gen_key_bucket(mesh.point(index));
            auto it = dPts.find(pt);
            if (it == dPts.end()) {
                faceIdx.push_back(dPts.size());

                tempPoly.push_back(dPts.size());

                dPts[pt] = dPts.size();
            } else {
                faceIdx.push_back(it->second);
                tempPoly.push_back(it->second);
            }
        }
        //- Check for problematic faces
        std::sort(faceIdx.begin(), faceIdx.end());
        auto it = std::unique(faceIdx.begin(), faceIdx.end());
        bool wasUnique = (it == faceIdx.end());

        if (wasUnique) {
            g["boundaries"].push_back({tempPoly});
        } else {
//            std::cerr << "Found duplicates!" << std::endl;
        }
    }
}