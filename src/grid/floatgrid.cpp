// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: floatgrid.cpp
// Implements floatgrid.hpp

#include "floatgrid.hpp"

using namespace fluidCore;
using namespace utilityCore;

floatgrid::floatgrid(){
	openvdb::initialize();
	vdbgrid = openvdb::FloatGrid::create(0.0f);
	type = VDB;
}

floatgrid::floatgrid(const gridtype& type, const glm::vec3& dimensions, const int& background){
	this->type = type;
	this->dimensions = dimensions;
	this->background = background;
	if(type==VDB){
		openvdb::initialize();
		vdbgrid = openvdb::FloatGrid::create(background);
	}else if(type==RAW){
		rawgrid = createGrid<float>(dimensions.x+1, dimensions.y+1, dimensions.z+1);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)dimensions.x+1),
			[=](const tbb::blocked_range<unsigned int>& r){
				for(unsigned int i=r.begin(); i!=r.end(); ++i){		
					for(unsigned int j=0; j<(int)dimensions.y+1; ++j){
						for(unsigned int k=0; k<(int)dimensions.z+1; ++k){
							rawgrid[i][j][k] = background;
						}
					}
				}
			}
		);
	}
}

floatgrid::floatgrid(openvdb::FloatGrid::Ptr newgrid){
	vdbgrid = newgrid;
	type = VDB;
}

floatgrid::~floatgrid(){
	if(type==RAW){
		deleteGrid<float>(rawgrid, dimensions.x+1, dimensions.y+1, dimensions.z+1);
	}else if(type==VDB){
		vdbgrid->clear();
		vdbgrid.reset();
	}
}

float floatgrid::getCell(const glm::vec3& index){
	return getCell((int)index.x, (int)index.y, (int)index.z);
}

float floatgrid::getCell(const int& x, const int& y, const int& z){
	float cell = background;
	if(type==VDB){
		openvdb::Coord coord = openvdb::Coord(x,y,z);
		openvdb::FloatGrid::Accessor accessor = vdbgrid->getAccessor();
		cell = accessor.getValue(coord);
	}else if(type==RAW){
		cell = rawgrid[x][y][z];
	}
	return cell;
}

void floatgrid::setCell(const glm::vec3& index, const float& value){
	setCell((int)index.x, (int)index.y, (int)index.z, value);
}

void floatgrid::setCell(const int& x, const int& y, const int& z, const float& value){
	if(type==VDB){
		SetCellLock.lock();
		{
			openvdb::Coord coord = openvdb::Coord(x,y,z);
			openvdb::FloatGrid::Accessor accessor = vdbgrid->getAccessor();
			accessor.setValue(coord, value);
		}
		SetCellLock.unlock();
	}else if(type==RAW){
		rawgrid[x][y][z] = value;
	}
}

float floatgrid::getInterpolatedCell(const glm::vec3& index){
	return getInterpolatedCell(index.x, index.y, index.z);
}

float floatgrid::getInterpolatedCell(const float& x, const float& y, const float& z){
	float value;
	GetInterpolatedCellLock.lock();
	{
		openvdb::Vec3f p(x,y,z);
		openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> interpolator(
													  vdbgrid->constTree(), vdbgrid->transform());
		value = interpolator.wsSample(p);
	}
	GetInterpolatedCellLock.unlock();
	return value;
}

openvdb::FloatGrid::Ptr& floatgrid::getVDBGrid(){
	return vdbgrid;
}

void floatgrid::writeVDBGridToFile(std::string filename){
	openvdb::io::File file(filename);
	openvdb::GridPtrVec grids;
    grids.push_back(vdbgrid);
    file.write(grids);
    file.close();
}

void floatgrid::clear(){
	if(type==VDB){
		vdbgrid->clear();
	}else if(type==RAW){
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)dimensions.x+1),
			[=](const tbb::blocked_range<unsigned int>& r){
				for(unsigned int i=r.begin(); i!=r.end(); ++i){
					for(unsigned int j=0; j<(int)dimensions.y+1; ++j){
						for(unsigned int k=0; k<(int)dimensions.z+1; ++k){
							rawgrid[i][j][k] = background;
						}
					}
				}
			}
		);	
	}
}
