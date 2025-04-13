# City4CFD

[![build](https://img.shields.io/github/actions/workflow/status/tudelft3d/City4CFD/build.yml?branch=main&style=flat-square)](https://github.com/tudelft3d/City4CFD/actions/workflows/build.yml)
[![docs](https://img.shields.io/badge/docs-Wiki-yellow?style=flat-square)](https://github.com/tudelft3d/City4CFD/wiki)
[![GitHub license](https://img.shields.io/github/license/tudelft3d/City4CFD?style=flat-square)](https://github.com/tudelft3d/City4CFD/blob/master/LICENSE)

![welcome_figure](/docs/images/workflow.png)

City4CFD--*City for CFD*--is a tool that aims to automatically reconstruct high-detailed 3D city geometries tailored for microscale urban flow simulations.

It can automatically create a terrain from a point cloud and imprint different surfaces (e.g. green areas, water, roads).

It enables us to reconstruct buildings from different sources and their combination, such as:
  - Reconstruction with the combination of 2D polygons and a point cloud,
  - Extrusion of footprints containing height or floor number attributes,
  - The import of existing building models.

The reconstruction using the combination of point clouds and 2D polygons can be done in [LoD2.2, LoD1.3, and LoD1.2](https://3d.bk.tudelft.nl/lod/). You can use the [complexity factor](https://github.com/tudelft3d/City4CFD/wiki/Features#buildings) to tune the overall complexity of the reconstructed geometry.

The resulting geometry is watertight -- buildings and surfaces are seamlessly integrated into a terrain.

It can automatically or manually define the zone of influence and domain boundaries.

If you happen to use it, feedback is very much appreciated.

The LoD2.2 and LoD1.3 reconstructions are based on [roofer](https://github.com/3DBAG/roofer). If you are interested in applications other than urban flow simulations, or you want to create your own 3D city modelling pipeline, we suggest checking out that project.

City4CFD is developed by the [3D Geoinformation Research Group](https://3d.bk.tudelft.nl/) at the Delft University of Technology.

## Data formats
**Point clouds** can be imported in LAS/LAZ, TXT/XYZ, or PLY format. We ask separately for ground and building points. While some datasets contain building-ground classification, some do not. Our [point cloud preparation tool](https://github.com/ipadjen/City4CFD_doc/wiki/Point-clouds#automatic-preparation) can extract ground and building points from user-defined classes, or use the [Cloth Simulation Filter](http://ramm.bnu.edu.cn/projects/CSF/) to separate the ground and non-ground points. If you would like to check your points, see if they are classified, or even conduct the filtering and classification yourself, we suggest you use [CloudCompare](https://www.danielgm.net/cc/).

**2D data** (polygons) are imported in [GDAL-supported formats](https://gdal.org/drivers/vector/index.html). For all pre-processing related to polygons you can use [QGIS](https://qgis.org/en/site/).

**Geometry import** supports the following formats: OBJ, STL, PLY, OFF, VTP, and CityJSON.

**Output** is in the following formats: OBJ, STL, and CityJSON. The ID of each polygon is preserved, and there is a 1-to-1 mapping between the input and the output.

## Installation
You can directly compile City4CFD on your system using cmake, run it through a Docker container, or install using Homebrew in the case of macOS.

### Build from source
The following libraries are required to build the project:
- [CGAL](https://www.cgal.org/) >= 6.0.1
- Boost >= 1.66
- Eigen >= 3.3.4
- GMP >= 4.2
- MPFR >= 2.2.1
- GDAL >= 3.0

*OpenMP* is an optional dependency.

GMP, MFPR, and Eigen are necessary dependencies for CGAL.

Dependencies are generally available in Linux distributions, e.g. in Debian/Ubuntu/Mint:
```
sudo apt-get install libmpfr-dev libgmp-dev libboost-all-dev libeigen3-dev libomp-dev libgdal-dev
```

CGAL can be directly downloaded from the [release page](https://github.com/CGAL/cgal/releases/tag/v6.0.1) (-library). No install is necessary, only the path to the unzipped folder is required (see below).

In macOS you can install all dependencies with Homebrew:

```
brew install cmake boost cgal eigen libomp gdal
```

The project uses CMake to generate makefiles, so make sure it is also installed.

To build City4CFD, do the following:
```
mkdir build && cd build
cmake .. -DCGAL_DIR=/path/to/cgal/dir
make
./city4cfd
```
In case of Homebrew, you do not have to use the ```-DCGAL_DIR``` parameter. You can speed up compilation by typing *make -j $numcores* where *$numcores* is the number of threads you can spare for compilation.

### Docker
We offer built [Docker](https://www.docker.com/) images for every release, available at the [Docker Hub](https://hub.docker.com/r/tudelft3d/city4cfd). Running [the docker script](https://github.com/tudelft3d/City4CFD/tree/main/docker/run) for the first time will pull the docker image from the Docker Hub repository.

### macOS
Mac users can install City4CFD through Homebrew:

```
brew install tudelft3d/software/city4cfd
```

## Getting started

The folder *examples* contains example datasets you can run for your first reconstruction. You can run your first reconstruction from the `/examples/TUD_Campus` folder by typing:
```
mkdir results
../../build/city4cfd config_bpg.json --output_dir results
```
in case of building from a source.

To run through a Docker container, you can use one of the scripts in ```docker/run/```. The script with the extension ```.sh``` can be used in Linux and macOS, the one with the extension ```.ps1``` in Windows Powershell, and the last one with ```.bat``` in Windows Command Prompt. You have to run  a script (you can copy it beforehand) from the root directory of the project (e.g. ```examples/TUD_Campus```), and the arguments are the same as for the compiled executable, e.g.: 

```
../../docker/run/city4cfd_run.sh city4cfd config_bpg.json --output_dir results
```

The script pulls the ```latest``` release from the Docker Hub. For a specific release, replace ```latest``` in the script with the released version tag, e.g. ```0.1.0```. In Linux systems, you will probably have to run the command as a sudo unless you create a 'docker' group and add users to it.

More information on the project can be found in the documentation.

## Documentation
The [wiki section](https://github.com/tudelft3d/City4CFD/wiki) of this project has details on reconstruction setup and also contains information and suggestions on data preparation.

## Citation
If you use City4CFD in a scientific context, please cite the following papers:

Ivan Pađen, Clara García-Sánchez, and Hugo Ledoux (2022). Towards Automatic Reconstruction of 3D City Models Tailored for Urban Flow Simulations. *Frontiers in Built Environment*, 8, 2022 [[DOI](https://doi.org/10.3389/fbuil.2022.899332)][[BibTeX](https://github.com/tudelft3d/City4CFD/blob/master/CITATION.bib)]

Ivan Pađen, Ravi Peters, Clara García-Sánchez, and Hugo Ledoux (2024). Automatic high-detailed building reconstruction workflow for urban microscale simulations. *Building and Environment*, 265, 2024 [[DOI](https://doi.org/10.1016/j.buildenv.2024.111978)][[BibTeX](https://github.com/tudelft3d/City4CFD/blob/master/CITATION.bib)]

## Acknowledgements
We would like to acknowledge the authors of the supporting libraries we use in this project:
[CGAL](https://github.com/CGAL/cgal), [CSF](https://github.com/jianboqi/CSF), [GDAL](https://github.com/OSGeo/gdal), [LAStools](https://github.com/LAStools), [nlohmann/json](https://github.com/nlohmann/json), [valijson](https://github.com/tristanpenman/valijson)
