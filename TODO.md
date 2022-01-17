# TODO
## Currently working
    - Import CityJSON buildings
        - Speed optimization

## High-priority
    - Averaged surfaces
        - Terrain smoothing must not touch averaged surface layers
        - Handling averaged surfaces that touch buildings
    - No terrain situation
    - Option to choose unique ID

## Longer-term
    - Reconstruct LoD1.2 from polygon height attribute
    - Manual extension of domain from building in m
    - Angled buffer zone
    - Polygon simplification
    - Option to translate the terrain up top
    - CityJSON overhaul

## Tested lib verisons
    - cmake 3.1
    - CGAL 5.0.2
    - Boost 1.65
