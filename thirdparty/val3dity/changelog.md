
# Changelog

## [Unreleased]
- validation of topological relationships between features, eg ensuring that buildings in a city do not overlap
- support for all GML3 primitives (for IndoorGML): the so-called "Compact Geometries" (http://schemas.opengis.net/gml/3.3/geometryCompact.xsd)


## [2.4.0] - 2023-06-27
### Added
- support for [CityJSON Lines (CityJSONL) files](https://www.cityjson.org/specs/#text-sequences-and-streaming-with-cityjsonfeature), the 1st line of the file must be with the `"transform"`, best is to use cjio output
- support for 3D types of [JSON-FG](https://github.com/opengeospatial/ogc-feat-geo-json), only the `Polyhedron` and `MultiPolyhedron` (not for `Prism` and `MultiPrism`) (JSON-FG is still beta, so I guess this might have to change)
- added a simple Python-based viewer to visualise the errors, see `/tools/viz3dity/` for more details
- API now supports OBJ and OFF; all formats supported by val3dity actually
### Changed
- the API is new and simpler, only 2 functions are now used: `is_valid()` + `validate()` 
- error 609 has been renamed 906 (no geometry is for more formats than only CityJSON)
- simplified the code for CDT of each surface (faster too)
- upgraded to latest version of pugixml, to parse XML files


## [2.3.1] - 2022-09-15
### Changed
- the structure of the report was fixed (it is more flat, as it should have always been), it made little sense for Composite/Multi/Solid
- geometry IDs in the report are 0-based (as they should, a bug made them start at 1, now fixed)
- compilation under Windows to get a binary is now automised with GitHub Actions


## [2.3.0] - 2022-03-29
### Added
- val3dity can now be compiled as a library and easily used in other C++ binary (it's GPLv3 license, watch out)
- this allows val3dity to be used with Python now, see https://github.com/tudelft3d/val3ditypy/
- unit test for IndoorGML files
- support for the simple format tu3djson (https://github.com/tudelft3d/tu3djson)
### Changed
- fixed a few bugs for the parsing of IndoorGML files
- LoDs of Buildings are now taken into account when validating the topological relationships between BuildingParts: the tests are now per LoD now (so no "clash" between LoD1.3 and LoD2.2 if stored in the same file)
### Removed
- removed the support for CityGML files and 'vanilla' GML files. Only IndoorGML files can be read, and support for GML 3.3 is not fully there yet. Most errors/issues that are reported with val3dity have to do with a faulty parser of CityGML, and not with val3dity itself. I cannot justify spending 100s of hours fixing the parser, I tried hard but this not possible anymore. If someone wants to do it then great I am willing to help them, but I will not invest time on this anymore. Just upgrade your files to CityJSON (https://www.cityjson.org/tutorials/conversion/), I'll always maintain a full CityJSON parser!


## [2.2.0] - 2020-05-14
### Added
- support for IndoorGML datasets, not just the geometries but specific validation tests for the primal-dual consistency were added
- support for GeometryTemplates in CityJSON (https://www.cityjson.org/specs/#geometry-templates)
### Changed
- completely revamped the JSON error report, it's better and clearer
- no more HTML output for the report (output of report is now with option `--report`, be careful), but you can just drag your JSON report to the HTML page (in folder `/report_browser/index.html`) or at http://geovalidation.bk.tudelft.nl/val3dity/browser/
- improved the error codes a bit: new errors for IndoorGML (7xx) + for handy faulty inputs (error 904)
- code must now be compiled against CGAL5, easier to install and compile for everyone
- CityGML support is deprecated. This is the last version that will support CityGML, I will remove it completely from v2.3+. Most errors/issues reported have to do with a faulty parser of CityGML, and not with val3dity itself. I cannot justify spending 100s of hours fixing the parser, I tried hard but this not possible anymore. If someone wants to do it then great I am willing to help them, but I will not invest time on this anymore. Just upgrade your files to CityJSON (https://www.cityjson.org/help/users/conversion/), I'll always maintain a full CityJSON parser!


## [2.1.1] - 2019-04-30
### Changed
- several small bug fixes: wrong handling of some GML namespaces is one of them
- tol_snap can't be disabled anymore: just too error-prone and means nothing. Give a super small value instead, like 1e-12


## [2.1.0] - 2018-02-14
### Changed
- the report is slightly changed: the CityObjects/Primitives subdivision in the report is replaced by Features having Primitives, where Features are the top-level objects, containing one or many Primitives. This is done to allow other data model to be used as input, eg IndoorGML is coming soon!
### Added
- pytest allows to specify the binary to use if not under 'build/'


## [2.0.4] - 2018-01-26
### Changed
- fixed bug where inner rings of gml:Polygon wasn't read properly when that polygon is referenced by an XLink.


## [2.0.3] - 2018-01-23
### Changed
- each surface are now triangulated with a constrained Delaunay triangulation (CDT) instead of a non-Delaunay one. That improves stability, especially when there are nearly collinear vertices in a surfaces: slivers are avoided, which means less errors 204, and less wrong errors caused by numerical stability.


## [2.0.2] - 2018-01-13 
### Changed
- fixes issue with wrong errors 104 often raised. Explanation is there: https://github.com/CGAL/cgal/issues/2733. This releases fixes it with a hack, but I think this is solid (touch wood).


## [2.0.1] - 2017-12-12 
### Changed
- fixed bug with least-square fitting that caused wrong 203 errors


## [2.0.0] - 2017-11-27 
### Added
- completely new reports: in JSON and HTML (browsable version)
- complete unit test suite (with pytest)
- support for all 3D primitives in ISO19107
- support for new formats: CityJSON, OFF


## [2.0 beta 1] - 2017-08-01 
### Added
- CompositeSolid, MultiSolid
- CityObjects in the report
- support for CityJSON
