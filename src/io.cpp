#include "io.h"

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

/*
void IO::output_cityjson(const TopoFeature* terrain, const std::vector<PolyFeature*>& buildings,
                         const TopoFeature* bnd) {
    std::ofstream of;
    nlohmann::json j;

    j["type"] = "CityJSON";
    j["version"] = "1.0";
    j["metadata"] = {};
    // bbox calc todo
//    double b[] = {bg::get<bg::min_corner, 0>(_bbox),
//                  bg::get<bg::min_corner, 1>(_bbox),
//                  0,
//                  bg::get<bg::max_corner, 0>(_bbox),
//                  bg::get<bg::max_corner, 1>(_bbox),
//                  0};
//    j["metadata"]["geographicalExtent"] = b;
    j["metadata"]["referenceSystem"] = "urn:ogc:def:crs:EPSG::7415";
    std::unordered_map< std::string, unsigned long > dPts;
    for (auto& f : buildings) {
        f->get_cityjson(j, dPts);
    }
    //-- vertices
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
    of << j.dump() << std::endl;

}
*/
