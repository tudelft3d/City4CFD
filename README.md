[![docs](https://img.shields.io/badge/docs-Wiki-brightgreen)](https://github.com/tudelft3d/City4CFD/wiki)
[![GitHub license](https://img.shields.io/github/license/tudelft3d/City4CFD)](https://github.com/tudelft3d/City4CFD/blob/master/LICENSE)

[//]: # ([![GitHub issues]&#40;https://img.shields.io/github/issues/tudelft3d/3dfier&#41;]&#40;https://github.com/tudelft3d/3dfier/issues&#41;)
[//]: # ([![DOI]&#40;https://joss.theoj.org/papers/10.21105/joss.02866/status.svg&#41;]&#40;https://doi.org/10.21105/joss.02866&#41;)


# City4CFD

City for CFD is a tool that aims to automatically reconstruct 3D city geometries tailored for microscale urban flow simulations.

It can create a terrain from a point cloud and imprint different surfaces (e.g. green areas, water, roads).

It enables the reconstruction of buildings from different sources and their combination, such as:
- Reconstruction with the combination of 2D polygons and a point cloud,
- Extrusion of footprints containing height or floor number attributes,
- The import of existing building models.

The resulting geometry is watertight -- buildings and surfaces are seamlessly integrated into a terrain.

It can automatically or manually define the zone of influence and domain boundaries.

If you happen to use it, feedback is very much appreciated.

City4CFD is developed by the [3D Geoinformation Research Group](https://3d.bk.tudelft.nl/) at the Delft University of Technology.

## Data formats
**Point clouds** can be imported in XYZ or PLY format. We ask separately for ground and building points. While some datasets contain building-ground classification, some do not. In case your LAS/LAZ file has buildings/ground classified, you can use [point cloud preparation script](https://github.com/tudelft3d/City4CFD/blob/main/tools/prepare_point_cloud.sh) to create separate files. If buildings and terrain are under the same class, or vegetation is not filtered out, we suggest you use [CloudCompare](https://www.danielgm.net/cc/) to prepare points.

**2D data** (polygons) are imported in GeoJSON format. For all pre-processing related to polygons, including conversion to GeoJSON, you can use [QGIS](https://qgis.org/en/site/).

**Geometry import** supports the following formats: OBJ, STL, PLY, OFF, VTP, and CityJSON.

**Output** is in the following formats: OBJ, STL, and CityJSON. The ID of each polygon is preserved, and there is a 1-to-1 mapping between the input and the output.

## Prerequisites

The following libraries are required to build the project:
- [CGAL](https://www.cgal.org/) version 5
- Boost >= 1.66

Both dependencies are generally available in Linux distributions (Debian/Ubuntu/Mint) as *libcgal-dev* and *libboost-dev*, and in macOS with Homebrew as *cgal* and *boost*.
The project uses CMake to generate makefiles, so make sure it is also installed.

## Installation

To build City4CFD, do the following:
```
mkdir build && cd build
cmake ..
make
./City4CFD
```
You can speed up compilation by typing *make -j $numcores* where *$numcores* is the number of threads you can spare for compilation.

## Getting started

The folder *examples* contains example datasets you can run for your first reconstruction. You can run your first reconstruction from the `/examples/TUD_Campus` folder by typing:
```
mkdir results
../../build/City4CFD config_bpg.json --output_dir results
```

More information on the project can be found in the documentation.

## Documentation
The [wiki section](https://github.com/tudelft3d/City4CFD/wiki) of this project has details on reconstruction setup and also contains information and suggestions on data preparation.


## Acknowledgements
We would like to acknowdledge the authors of supporting libraries we use in this project:
[CGAL](https://github.com/CGAL/cgal), [nlohmann/json](https://github.com/nlohmann/json), [valijson](https://github.com/tristanpenman/valijson), [LAStools](https://github.com/LAStools)
