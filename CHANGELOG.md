# Changelog
## [0.6.3] 2025-03-07
### Fixed
- CityJSON terrain output

## [0.6.2] - 2025-02-27
### Changed
- Improved inserting surface layer polygons when terrain pc is missing
### Fixed
- Minor bugfixes

## [0.6.1] - 2025-02-06
### Fixed
- Docker image hotfix
- CGAL version update in CMakeLists

## [0.6.0] - 2025-01-03
### Changed
- Updated to CGAL 6

## [0.5.0] - 2024-07-29 - breaking changes
### Added
- LoD2.2 and LoD1.3 reconstruction
- Multiple reconstruction regions
- Geometric validity check
- Fallbacks for failed reconstructions
- Polygon simplification tool

## [0.4.6] - 2024-03-11
### Changed
- Improved robustness in handling bad polygons
### Fixed
- Code crashing when flattening specific cases

## [0.4.5] - 2024-02-12
### Fixed
- Dependency hotfix

## [0.4.4] - 2024-01-13
### Changed
- Improved fallbacks in case building import fails in late stages
- Improved polygon flattening when buildings are adjacent
### Fixed
- Minor bugfixes

## [0.4.3] - 2023-08-25
### Fixed
- Issue with GDAL on Ubuntu 20.04

## [0.4.2]  - 2023-07-12
### Fixed
- Problematic compilation in debug mode

## [0.4.1]  - 2023-06-05
### Fixed
- STL export bugfix

## [0.4.0] - 2023-05-03
### Added
- (breaking) Support for many polygon formats (through GDAL)
- Optional flag that enforces buildings and terrain intersection
### Changed
- Updated LAStools to v2.0.2
- Minor bugfixes

## [0.3.0] - 2023-01-18
### Added
- Building surface refinement
- Sharper boundaries of flattened polygons
- Minimum building height as an argument
### Fixed
- Bad reconstruction of buildings with low height
- Minor bugfixes

## [0.2.0] - 2022-10-07
### Added
- Point cloud preparation tool city4cfd_pcprep
### Changed
- (breaking) Terrain smoothing overhaul
### Removed
- Previous point cloud preparation script

## [0.1.2] - 2022-09-13
### Changed
- Bad triangles handling

## [0.1.1] - 2022-09-02
### Fixed
- Bad triangles at surface edges
- Polygons close to the boundary are defined as 'out of bounds'

## [0.1.0] - 2022-08-14 
First release
