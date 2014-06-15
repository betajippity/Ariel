// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: intgrid.cpp
// Implements intgrid.hpp

#include "intgrid.hpp"

namespace fluidCore{

intgrid::intgrid(const gridtype& type, const glm::vec3& dimensions, const int& background){
	this->type = type;
	this->dimensions = dimensions;
	this->background = background;
	if(type==VDB){
		openvdb::initialize();
		vdbgrid = openvdb::Int32Grid::create(background);
	}else if(type==RAW){
		rawgrid = createGrid<int>(dimensions.x+1, dimensions.y+1, dimensions.z+1);

		tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)dimensions.x+1),
			[=](const tbb::blocked_range<unsigned int>& r){
				for(unsigned int i=r.begin(); i!=r.end(); ++i){	
					for(int j=0; j<(int)dimensions.y+1; ++j){
						for(int k=0; k<(int)dimensions.z+1; ++k){
							rawgrid[i][j][k] = background;
						}
					}
				}
			}
		);
	}
}

intgrid::~intgrid(){
	if(type==RAW){
		deleteGrid<int>(rawgrid, dimensions.x+1, dimensions.y+1, dimensions.z+1);
	}else if(type==VDB){
		vdbgrid->clear();
		vdbgrid.reset();
	}
}

int intgrid::getCell(const glm::vec3& index){
	return getCell((int)index.x, (int)index.y, (int)index.z);
}

int intgrid::getCell(const int& x, const int& y, const int& z){
	int cell = background;
	if(type==VDB){
		openvdb::Coord coord = openvdb::Coord(x,y,z);
		openvdb::Int32Grid::Accessor accessor = vdbgrid->getAccessor();
		cell = accessor.getValue(coord);
	}else if(type==RAW){
		cell = rawgrid[x][y][z];
	}
	return cell;
}

void intgrid::setCell(const glm::vec3& index, const int& value){
	setCell((int)index.x, (int)index.y, (int)index.z, value);
}

void intgrid::setCell(const int& x, const int& y, const int& z, const int& value){
	if(type==VDB){
		SetCellLock.lock();
		{
			openvdb::Coord coord = openvdb::Coord(x,y,z);
			openvdb::Int32Grid::Accessor accessor = vdbgrid->getAccessor();
			accessor.setValue(coord, value);
		}
		SetCellLock.unlock();
	}else if(type==RAW){
		rawgrid[x][y][z] = value;
	}
}

openvdb::Int32Grid::Ptr& intgrid::getVDBGrid(){
	return vdbgrid;
}

void intgrid::writeVDBGridToFile(std::string filename){
	openvdb::io::File file(filename);
	openvdb::GridPtrVec grids;
    grids.push_back(vdbgrid);
    file.write(grids);
    file.close();
}

gridtype intgrid::getGridType(){
	return type;
}
}
