#ifndef CITYCFD_IO_H
#define CITYCFD_IO_H

#include "iostream"
#include "definitions.h"

// TODO; hardcoded for now
struct ConfigData {
    Point_2     pointOfInterest = Point_2(85420,446221);
    double radiusOfInterest = 350.0;
    double dimOfDomain = 1000.0;
    double topHeight   = 100.0;
    // note: handle when radiusOfInterst is larger than dimOfDomain
};


void          read_config();
std::string   gen_key_bucket(const Point_3* p);
void          get_obj(const Mesh& mesh, std::string& fs, std::string& bs);
void          get_obj(const Mesh& mesh, std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts);
void          get_obj(const Mesh& mesh, std::string& fs, std::string& bs, std::string& className);

#endif //CITYCFD_IO_H
