#include "io.h"

void read_config() {

}

std::string gen_key_bucket(const Point_3* p) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << p->x() << " " << p->y() << " " << std::setprecision(2) << p->z();
    return ss.str();
}

//-- Output functions
void output_obj(const TopoFeature* terrain, const std::vector<PolyFeature*>& buildings, const TopoFeature* boundary, bool outputSeparately) {
    std::unordered_map<std::string, unsigned long> dPts;
    std::vector<std::ofstream> of;
    std::vector<std::string>   fs, bs;

    int count = 0;
    if (outputSeparately) {
        of.reserve(4);
        of.emplace_back().open("Terrain.obj"); fs.emplace_back(); bs.emplace_back();
    } else {
        of.emplace_back().open("Mesh.obj"); fs.emplace_back(); bs.emplace_back();
    }
    //-- Output terrain
    bs[count] += "\ng Terrain";
    get_obj(terrain->get_mesh(), fs[count], bs[count], dPts);

    if (outputSeparately) {
        ++count; dPts.clear();
        of.emplace_back().open("Building.obj"); fs.emplace_back(); bs.emplace_back();
    }
    //-- Output buildings
    bs[count] += "\ng Building";
    for (auto& f : buildings) {
        if (!f->is_active()) continue;
        get_obj(f->get_mesh(), fs[count], bs[count], dPts);
    }

    if (outputSeparately) {
        ++count; dPts.clear();
        of.emplace_back().open("Sides.obj"); fs.emplace_back(); bs.emplace_back();
    }
    //-- Output side and top boundary
    bs[count] += "\ng Sides";
    get_obj(boundary->get_mesh(), fs[count], bs[count], dPts);

    if (outputSeparately) {
        ++count; dPts.clear();
        of.emplace_back().open("Top.obj"); fs.emplace_back(); bs.emplace_back();
    }
    bs[count] += "\ng Top";
    get_obj(dynamic_cast<const Boundary*>(boundary)->get_top_mesh(), fs[count], bs[count], dPts);

    //-- Write to file
    int i = 0;
    do {
        of[i] << fs[i] << bs[i];
        of[i].close();
    } while (i++ < count);
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