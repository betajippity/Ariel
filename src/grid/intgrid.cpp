// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: intgrid.cpp
// Implements intgrid.cpp

#include "intgrid.hpp"

using namespace fluidCore;
using namespace utilityCore;

intgrid::intgrid(const int& background){
	openvdb::initialize();
	grid = openvdb::Int32Grid::create(background);
}

intgrid::~intgrid(){

}

int intgrid::getCell(const vec3& index){
	return getCell((int)index.x, (int)index.y, (int)index.z);
}

int intgrid::getCell(const int& x, const int& y, const int& z){
	openvdb::Coord coord = openvdb::Coord(x,y,z);

	openvdb::Int32Grid::Accessor accessor = grid->getAccessor();
	return accessor.getValue(coord);
}

void intgrid::setCell(const vec3& index, const int& value){
	setCell((int)index.x, (int)index.y, (int)index.z, value);
}

void intgrid::setCell(const int& x, const int& y, const int& z, const int& value){
	openvdb::Coord coord = openvdb::Coord(x,y,z);

	openvdb::Int32Grid::Accessor accessor = grid->getAccessor();

	if(epsilonCheck(length(value), 0.0f)==true){
		accessor.setValueOff(coord, value);
	}else{
		accessor.setValue(coord, value);
	}
}
