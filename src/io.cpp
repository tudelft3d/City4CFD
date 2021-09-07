#include "io.h"

void read_config() {

}

std::string gen_key_bucket(const Point_3* p) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << p->x() << " " << p->y() << " " << std::setprecision(2) << p->z();
    return ss.str();
}

void get_obj(const Mesh& mesh, std::string& fs, std::string& bs) {
    for (auto &vert : mesh.vertices()) {
        fs += "\nv " + gen_key_bucket(&mesh.point(vert));
    }

    for (auto& face : mesh.faces()) {
        bs += "\nf";
        for (auto index : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            bs += " " + std::to_string(index.idx() + 1);
        }
    }

}

void get_obj(const Mesh& mesh,
             std::string& fs,
             std::string& bs,
             std::unordered_map<std::string, unsigned long>& dPts)
{
    for (auto& face : mesh.faces()) {
        std::vector<unsigned long> faceIdx;
        std::string fsTemp;
        std::string bsTemp;
        for (auto index : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            std::string pt = gen_key_bucket(&mesh.point(index));
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
        sort(faceIdx.begin(), faceIdx.end());
        auto it = std::unique( faceIdx.begin(), faceIdx.end() );
        bool wasUnique = (it == faceIdx.end() );

        if (wasUnique) {
            bs += "\nf";
            bs += bsTemp;
        } else {
//            std::cerr << "Found duplicates!" << std::endl;
        }
    }
}

void get_obj(const Mesh &mesh, std::string &fs, std::string &bs, std::string &className) {
    fs += "o " + className;
    for (auto &vert : mesh.vertices()) {
        fs += "\nv " + gen_key_bucket(&mesh.point(vert));
    }

    for (auto &face : mesh.faces()) {
        bs += "\nf";
        for (auto index : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            bs += " " + std::to_string(index.idx() + 1);
        }
    }
}
