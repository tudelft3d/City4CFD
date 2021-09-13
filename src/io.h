#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "definitions.h"
#include "config.h"

namespace IO {
    //-- Input functions

    //-- Output functions
    void get_obj_pts(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);
    void get_stl_pts(Mesh& mesh, std::string& fs);

//    void output_cityjson(const TopoFeature* terrain, const std::vector<PolyFeature*> &buildings, const TopoFeature* bnd);

    //-- Templated function
    template <typename T>
    std::string gen_key_bucket(const T& p) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << p.x() << " " << p.y() << " " << p.z();
        return ss.str();
    }

}

#endif //CITYCFD_IO_H