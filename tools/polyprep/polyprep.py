"""
  polyprep

  Copyright (c) 2024, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
"""

import argparse
from shapely.geometry import shape, mapping, MultiPolygon, Polygon
from shapely.ops import unary_union
from shapely.prepared import prep
import fiona
from rtree import index

def close_holes_in_polygons(polygons):
    """
    Removes holes from the provided list of polygons.
    """
    def close_holes(polygon):
        if isinstance(polygon, Polygon):
            return Polygon(polygon.exterior)
        elif isinstance(polygon, MultiPolygon):
            return MultiPolygon([Polygon(p.exterior) for p in polygon.geoms])
        else:
            return polygon

    return [close_holes(polygon) for polygon in polygons]

def determine_properties_from_components(components, polygon_data):
    """
    Transfers the properties of the polygon with the largest area in case of merging polygons.
    """
    properties_values = []

    idx = index.Index()
    for i, (geom, props) in enumerate(polygon_data):
        idx.insert(i, geom.bounds)

    if isinstance(components, MultiPolygon):
        components = list(components.geoms)
    else:
        components = [components]

    for component in components:
        overlapping_polygons = []

        for i in idx.intersection(component.bounds):
            data_geom, data_props = polygon_data[i]

            prep_data_geom = prep(data_geom)

            if prep_data_geom.intersects(component):
                overlapping_polygons.append((data_geom, data_props))

        if not overlapping_polygons:
            properties_values.append(None)
            continue

        largest_polygon = max(overlapping_polygons, key=lambda x: x[0].area)
        properties_values.append(largest_polygon[1])

    return properties_values

def process_polygons(input_filename, output_filename, buffer_size, apply_convex_hull=False, remove_holes=0,
                     simplification_tol=0.):
    with fiona.open(input_filename, 'r') as input_src:
        input_schema = input_src.schema.copy()
        input_crs = input_src.crs

        output_schema = {
            'geometry': 'Polygon',
            'properties': input_schema['properties']
        }

        polygon_data = []
        for feature in input_src:
            geom = shape(feature['geometry'])
            properties = feature['properties']

            if isinstance(geom, Polygon):
                polygon_data.append((geom, properties))
            elif isinstance(geom, MultiPolygon):
                for poly in geom.geoms:
                    polygon_data.append((poly, properties))

        print("Buffering polygons...")
        buffered_polygons = [data[0].buffer(buffer_size, cap_style=3, join_style=2) for data in polygon_data]
        if remove_holes != 0:
            buffered_polygons = close_holes_in_polygons(buffered_polygons)

        print("Dissolving polygons...")
        dissolved_polygons = unary_union(buffered_polygons).buffer(-buffer_size, buffer_size, cap_style=3, join_style=2)

        if remove_holes == 2:  # additional pass to remove any new holes created
            if isinstance(dissolved_polygons, MultiPolygon):
                dissolved_polygons = MultiPolygon(close_holes_in_polygons(dissolved_polygons.geoms))
            else:
                dissolved_polygons = close_holes_in_polygons([dissolved_polygons])[0]

        new_properties_values = determine_properties_from_components(dissolved_polygons, polygon_data)

        polygons_to_write = []

        if apply_convex_hull:
            print("Extracting convex hull...")
            if isinstance(dissolved_polygons, MultiPolygon):
                convex_hulls = [polygon.convex_hull for polygon in dissolved_polygons.geoms]
            else:
                convex_hulls = [dissolved_polygons.convex_hull]

            polygon_data = list(zip(convex_hulls, new_properties_values))
            merged_hulls = unary_union(convex_hulls) # chulls need an extra dissolve pass as it can introduce overlaps
            new_properties_values = determine_properties_from_components(merged_hulls, polygon_data)

            if isinstance(merged_hulls, MultiPolygon):
                polygons_to_write = list(zip(list(merged_hulls.geoms), new_properties_values))
            else:
                polygons_to_write = list(zip([merged_hulls], new_properties_values))
        else:
            if isinstance(dissolved_polygons, MultiPolygon):
                polygons_to_write = list(zip(list(dissolved_polygons.geoms), new_properties_values))
            else:
                polygons_to_write = list(zip([dissolved_polygons], new_properties_values))

        print("Writing new polygons to file...")
        with fiona.open(output_filename, 'w', 'GeoJSON', output_schema, crs=input_crs) as output_src:
            for polygon, properties in polygons_to_write:
                final_properties = dict(input_src[0]['properties'])
                final_properties.update(properties)

                if simplification_tol > 0.:
                    polygon = polygon.simplify(simplification_tol)

                feature = {'geometry': mapping(polygon), 'properties': final_properties}
                output_src.write(feature)

        print("End")

def main():
    parser = argparse.ArgumentParser(
        description="polyprep: City4CFD polygon generalisation tool\n\n"
                    "This tool generalises GeoJSON polygon file so that subsequent building reconstruction results in "
                    "simplified geometries.\n"
                    "It can buffer, simplify, and optionally remove holes from polygons, among other features.\n\n"
                    "Example usage:\n"
                    "  python polyprep.py input.geojson output.geojson 1.0 --remove_holes 1 --simplification_tol 0.1",
        epilog="This program is released under the GNU Affero General Purpose Licence (AGPL) v3\n"
               "Copyright(c) 2024, 3D Geoinformation Research Group, TU Delft",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("input_filename", help="input filename (GeoJSON format)")
    parser.add_argument("output_filename", help="output filename (GeoJSON format)")
    parser.add_argument("buffer_size", type=float, help="buffer size for polygon merging")
    parser.add_argument("--convex_hull", action="store_true", help="simplify polygons by applying convex hull")
    parser.add_argument("--remove_holes", type=int, choices=[0, 1, 2], default=0, help="remove holes from polygons: 0 - No, 1 - Yes, 2 - Yes with additional pass after poligon merging")
    parser.add_argument("--simplification_tol", type=float, default=0.0, help="tolerance for polygon simplification using Douglas-Pueckler algrithm")
    parser.add_argument("--version", action="version", version="polyprep version 0.1.0")

    args = parser.parse_args()

    process_polygons(
        input_filename=args.input_filename,
        output_filename=args.output_filename,
        buffer_size=args.buffer_size,
        apply_convex_hull=args.convex_hull,
        remove_holes=args.remove_holes,
        simplification_tol=args.simplification_tol
    )

if __name__ == "__main__":
    main()
