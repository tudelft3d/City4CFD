#include "io.h"

void IO::output_obj(std::vector<TopoFeature*>& allFeatures) {
    using namespace config;
    std::vector<std::ofstream> of;
    std::vector<std::string>   fs(topoClassName.size()), bs(topoClassName.size());

    std::vector<std::unordered_map<std::string, unsigned long>> dPts(topoClassName.size());
    //-- Output points
    for (auto& f : allFeatures) {
        if (!f->is_active()) continue;
        if (outputSeparately)
            IO::get_obj_pts(f->get_mesh(), fs[f->get_class()], bs[f->get_class()], dPts[f->get_class()]);
        else
            IO::get_obj_pts(f->get_mesh(), fs[f->get_class()], bs[f->get_class()], dPts[0]);
    }

    //-- Add class name and output to file
    if (!outputSeparately) of.emplace_back().open(fileName + ".obj");
    for (int i = 0; i < fs.size(); ++i) {
        if (bs[i].empty()) continue;
        if (outputSeparately) of.emplace_back().open(fileName + "_" + topoClassName.at(i) + ".obj");

        of.back() << fs[i] << "\ng " << topoClassName.at(i) << bs[i];
    }
    for (auto& f : of) f.close();
}

void IO::output_stl(std::vector<TopoFeature*>& allFeatures) {
    using namespace config;
    std::vector<std::ofstream> of;
    std::vector<std::string>   fs(topoClassName.size());

    //-- Get all triangles
    for (auto& f : allFeatures) {
        if (!f->is_active()) continue;
        IO::get_stl_pts(f->get_mesh(), fs[f->get_class()]);
    }

    //-- Add class name and output to file
    if (!outputSeparately) of.emplace_back().open(fileName + ".stl");
    for (int i = 0; i < fs.size(); ++i) {
        if (fs[i].empty()) continue;
        if (outputSeparately) of.emplace_back().open(fileName + "_" + topoClassName.at(i) + ".stl");

        of.back() << "\nsolid " << topoClassName.at(i);
        of.back() << fs[i];
        of.back() << "\nendsolid " << topoClassName.at(i);
    }
    for (auto& f : of) f.close();
}

void IO::output_cityjson(std::vector<TopoFeature*>& allFeatures) {
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
        // TESTING ONLY BUILDINGS NOW
        if (f->get_class() != BUILDING) continue;
        //-- Get feature info
        nlohmann::json b;
        f->get_cityjson_info(b);

        //-- Get feature geometry
        nlohmann::json g;
        IO::get_cityjson_geom(f->get_mesh(), g, dPts, f->get_cityjson_primitive());

        //-- Append to main json struct
        b["geometry"].push_back(g);
        j["CityObjects"][f->get_id()] = b;
    }

    //-- Vertices
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

    of.open(fileName + ".json");
    of << j.dump() << std::endl;
}

//-- Output functions
void IO::get_obj_pts(const Mesh& mesh,
                     std::string& fs,
                     std::string& bs,
                     std::unordered_map<std::string, unsigned long>& dPts)
{
    for (auto& face : mesh.faces()) {
        std::vector<unsigned long> faceIdx;
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
        std::vector<unsigned long> faceIdx;
        std::vector<unsigned long> tempPoly;
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