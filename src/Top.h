#ifndef CITY4CFD_TOP_H
#define CITY4CFD_TOP_H

#include "Boundary.h"

class Top : public Boundary {
public:
    Top(const int outputLayerID);
    ~Top();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

#endif //CITY4CFD_TOP_H