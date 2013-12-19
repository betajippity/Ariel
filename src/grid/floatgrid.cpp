// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: floatgrid.cpp
// Implements floatgrid.hpp

#include "floatgrid.hpp"

using namespace fluidCore;
using namespace utilityCore;

floatgrid::floatgrid(){
	openvdb::initialize();
	grid = openvdb::FloatGrid::create(0.0f);
}

floatgrid::floatgrid(const float& background){
	openvdb::initialize();
	grid = openvdb::FloatGrid::create(background);
}

floatgrid::floatgrid(openvdb::FloatGrid::Ptr newgrid){
	grid = newgrid;
}

floatgrid::~floatgrid(){
	grid->clear();
	grid.reset();
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
	#pragma omp critical
	{
		openvdb::Coord coord = openvdb::Coord(x,y,z);

		openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

		accessor.setValue(coord, value);
	}
}

float floatgrid::getInterpolatedCell(const vec3& index){
	return getInterpolatedCell(index.x, index.y, index.z);
}

float floatgrid::getInterpolatedCell(const float& x, const float& y, const float& z){
	openvdb::Vec3f p(x,y,z);
	openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> interpolator(
														 grid->constTree(), grid->transform());
	return interpolator.wsSample(p);
}

openvdb::FloatGrid::Ptr& floatgrid::getVDBGrid(){
	return grid;
}

void floatgrid::writeVDBGridToFile(string filename){
	openvdb::io::File file(filename);
	openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();
}

void floatgrid::clear(){
	grid->clear();
}
