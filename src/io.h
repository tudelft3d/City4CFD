#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "definitions.h"
#include "TopoFeature.h"
#include "Boundary.h"

//-- Input functions


//-- Output functions
void output_obj(const TopoFeature* terrain, const std::vector<PolyFeature*>& buildings, const TopoFeature* bnd, bool outputSeparately);

void          read_config();
std::string   gen_key_bucket(const Point_3* p);
void          get_obj(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);

#endif //CITYCFD_IO_H