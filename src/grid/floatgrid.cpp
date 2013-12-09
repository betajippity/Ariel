// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: floatgrid.cpp
// Implements floatgrid.cpp

#include "floatgrid.hpp"

using namespace fluidCore;
using namespace utilityCore;

floatgrid::floatgrid(const float& background){
	openvdb::initialize();
	grid = openvdb::FloatGrid::create(background);
}

floatgrid::~floatgrid(){

}

float floatgrid::getCell(const vec3& index){
	return getCell((int)index.x, (int)index.y, (int)index.z);
}

float floatgrid::getCell(const int& x, const int& y, const int& z){
	openvdb::Coord coord = openvdb::Coord(x,y,z);

	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	return accessor.getValue(coord);
}

void floatgrid::setCell(const vec3& index, const float& value){
	setCell((int)index.x, (int)index.y, (int)index.z, value);
}

void floatgrid::setCell(const int& x, const int& y, const int& z, const float& value){
	openvdb::Coord coord = openvdb::Coord(x,y,z);

	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

	if(epsilonCheck(length(value), 0.0f)==true){
		accessor.setValueOff(coord, value);
	}else{
		accessor.setValue(coord, value);
	}
}
