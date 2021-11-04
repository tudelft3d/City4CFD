#ifndef CITYCFD_SIDES_H
#define CITYCFD_SIDES_H

#include "Boundary.h"

class Sides : public Boundary {
public:
    Sides(const int outputLayerID);
    ~Sides();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

#endif //CITYCFD_SIDES_H
