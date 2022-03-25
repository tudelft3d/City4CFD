#ifndef CITY4CFD_SIDES_H
#define CITY4CFD_SIDES_H

#include "Boundary.h"

class Sides : public Boundary {
public:
    Sides(const int outputLayerID);
    ~Sides();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

#endif //CITY4CFD_SIDES_H
