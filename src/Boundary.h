#ifndef CITYCFD_BOUNDARY_H
#define CITYCFD_BOUNDARY_H

#include "definitions.h"
#include "geomtools.h"
#include "io.h"
#include "TopoFeature.h"

class Boundary : public TopoFeature {
public:
    Boundary() = default;
    //TODO: get domain info from the config file
    Boundary(const ConfigData& configData);
    ~Boundary() = default;

    void set_bounds_to_pc(Point_set_3& pointCloud) const;
    void set_bounds_to_cdt(CDT& cdt) const;
    void add_buffer(Point_set_3& pointCloud);
    void get_mesh();

    TopoClass get_class() const override;
    void      output_feature(std::string& fs, std::string& bs, std::unordered_map<std::string, unsigned long>& dPts) const override;

private:
    //-- HARDCODED data, should be there from the config--//
    const Point_2 pointOfInterest        = Point_2(85420,446221);
    const double  radiusOfInterestConfig = 350.;
    const double  dimOfDomainConfig      = 1000.;
    const double  topHeightConfig        = 200.;
    //------------------//

    double               _dimOfDomain = dimOfDomainConfig;
    std::vector<Point_3> _outerPts;
    Mesh                 _meshTop; // _mesh from TopoFeature is for the sides
    double               _topHeight = topHeightConfig;
};


#endif //CITYCFD_BOUNDARY_H
