#include "io.h"

//-- Output functions
void IO::output_obj(const TopoFeature* terrain, const std::vector<PolyFeature*>& buildings, const TopoFeature* bnd) {
    using namespace config;
    std::unordered_map<std::string, unsigned long> dPts;
    std::vector<std::ofstream> of;
    std::vector<std::string>   fs, bs;

    //-- Output terrain
    fs.emplace_back(); bs.emplace_back();
    bs.back() += "\ng Terrain";
    get_obj(terrain->get_mesh(), fs.back(), bs.back(), dPts);

    if (outputSeparately) dPts.clear();
    //-- Output buildings
    fs.emplace_back(); bs.emplace_back();
    bs.back() += "\ng Building";
    for (auto& f : buildings) {
        if (!f->is_active() || f->get_class() != BUILDING) continue;
        get_obj(f->get_mesh(), fs.back(), bs.back(), dPts);
    }

    if (outputSeparately) dPts.clear();
    //-- Output side and top boundary
    fs.emplace_back(); bs.emplace_back();
    bs.back() += "\ng Sides";
    get_obj(bnd->get_mesh(), fs.back(), bs.back(), dPts);

    if (outputSeparately) dPts.clear();
    fs.emplace_back(); bs.emplace_back();
    bs.back() += "\ng Top";
    get_obj(dynamic_cast<const Boundary*>(bnd)->get_top_mesh(), fs.back(), bs.back(), dPts);

    //-- Write to file
    if (config::outputSeparately) {
        of.emplace_back().open("Terrain.obj");   of.back() << fs[0] << bs[0];
        of.emplace_back().open("Buildings.obj"); of.back() << fs[1] << bs[1];
        of.emplace_back().open("Sides.obj");     of.back() << fs[2] << bs[2];
        of.emplace_back().open("Top.obj");       of.back() << fs[3] << bs[3];
    } else {
        of.emplace_back().open("Mesh.obj");
        for (int i = 0; i < fs.size(); ++i) {
            of.back() << fs[i] << bs[i];
        }
    }
    for (auto& f: of) {
        f.close();
    }
}

void IO::get_obj(const Mesh& mesh,
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

void IO::output_stl(TopoFeature* terrain, std::vector<PolyFeature*>& buildings, TopoFeature* bnd) {
    std::vector<std::ofstream> of;
    std::vector<std::string>   fs;

    //-- Output terrain
    fs.emplace_back();
    fs.back() += "solid Terrain";
    get_stl(terrain->get_mesh(), fs.back());
    fs.back() += "\nendsolid Terrain";

    //-- Output buildings
    fs.emplace_back();
    fs.back() += "solid Building";
    for (auto& f :buildings) {
        if (!f->is_active() || f->get_class() != BUILDING) continue;
        get_stl(f->get_mesh(), fs.back());
    }
    fs.back() += "\nendsolid Building";

    //-- Output domain edges
    fs.emplace_back();
    fs.back() += "solid Sides";
    get_stl(bnd->get_mesh(), fs.back());
    fs.back() += "\nendsolid Sides";

    fs.emplace_back();
    fs.back() += "solid Top";
    get_stl(dynamic_cast<Boundary*>(bnd)->get_top_mesh(), fs.back());
    fs.back() += "\nendsolid Top";

    //-- Write to file
    if (config::outputSeparately) {
        of.emplace_back().open("Terrain.stl");   of.back() << fs[0];
        of.emplace_back().open("Buildings.stl"); of.back() << fs[1];
        of.emplace_back().open("Sides.stl");     of.back() << fs[2];
        of.emplace_back().open("Top.stl");       of.back() << fs[3];
    } else {
        of.emplace_back().open("Mesh.stl");
        for (auto& f : fs) {
            of.back() << f << "\n";
        }
    }
    for (auto& f: of) {
        f.close();
    }

}

void IO::get_stl(Mesh& mesh, std::string& fs) {
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