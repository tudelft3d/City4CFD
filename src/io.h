#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "definitions.h"
#include "config.h"
#include "TopoFeature.h"
#include "Boundary.h"

namespace IO {
//-- Input functions
    void read_config();


//-- Output functions
    void output_obj(const TopoFeature* terrain, const std::vector<PolyFeature*> &buildings, const TopoFeature* bnd);
    std::string gen_key_bucket(const Point_3* p);
    void get_obj(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);
}

#endif //CITYCFD_IO_H