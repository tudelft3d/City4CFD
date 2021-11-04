#ifndef CITYCFD_TOP_H
#define CITYCFD_TOP_H

#include "Boundary.h"

class Top : public Boundary {
public:
    Top(const int outputLayerID);
    ~Top();

    virtual void reconstruct() override;

    virtual TopoClass   get_class() const override;
    virtual std::string get_class_name() const override;
};

#endif //CITYCFD_TOP_H