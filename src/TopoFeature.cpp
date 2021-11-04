#include "TopoFeature.h"

//-- TopoFeature class
TopoFeature::TopoFeature()
    : _mesh(), _id(), _f_active(true), _outputLayerID(-1) {}

TopoFeature::TopoFeature(std::string pid)
    : _mesh(), _id(std::move(pid)), _f_active(true), _outputLayerID(-1) {}

TopoFeature::TopoFeature(int outputLayerID)
    : _mesh(), _id(), _f_active(true), _outputLayerID(outputLayerID) {
    if (_outputLayerID  >= _numOfOutputLayers) _numOfOutputLayers = _outputLayerID + 1;
}

TopoFeature::~TopoFeature() = default;

int TopoFeature::_numOfOutputLayers = 0;

int TopoFeature::get_num_output_layers() {
    return _numOfOutputLayers;
}

Mesh& TopoFeature::get_mesh() {
    return _mesh;
}

const Mesh& TopoFeature::get_mesh() const {
    return _mesh;
}

void TopoFeature::set_id(unsigned long id) {
    _id = std::to_string(id);
}

std::string TopoFeature::get_id() const {
    return _id;
}

const int TopoFeature::get_output_layer_id() const {
    return _outputLayerID;
}

bool TopoFeature::is_active() const {
    return _f_active;
}

void TopoFeature::deactivate() {
    _f_active = false;
}

void TopoFeature::get_cityjson_info(nlohmann::json& j) const {
    //TEMP UNTIL ALL FUNCTIONS ARE IMPLEMENTED
}

void TopoFeature::get_cityjson_semantics(nlohmann::json& g) const {
    // TEMP until I figure what to do with this
}


std::string TopoFeature::get_cityjson_primitive() const {
    //TEMP UNTIL ALL FUNCTIONS ARE IMPLEMENTED
    return "Nope";
}