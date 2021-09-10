#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "definitions.h"
#include "config.h"
#include "TopoFeature.h"
#include "Boundary.h"

namespace IO {
    //-- Input functions

    //-- Output functions
    void output_obj(const TopoFeature* terrain, const std::vector<PolyFeature*> &buildings, const TopoFeature* bnd);
    void get_obj(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);
    void output_stl(TopoFeature* terrain, std::vector<PolyFeature*> &buildings, TopoFeature* bnd);
    void get_stl(Mesh& mesh, std::string& fs);

    //-- Templated function
    template <typename T>
    std::string gen_key_bucket(const T& p) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << p.x() << " " << p.y() << " " << p.z();
        return ss.str();
    }

}

#endif //CITYCFD_IO_H