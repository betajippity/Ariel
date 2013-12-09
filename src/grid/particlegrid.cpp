// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegrid.cpp
// Implements particlegrid.hpp

#include "particlegrid.hpp"

using namespace fluidCore;

particlegrid::particlegrid(const vec3& dim){
	init((int)dim.x, (int)dim.y, (int)dim.z);
}

particlegrid::particlegrid(const int& x, const int& y, const int& z){
	init(x, y, z);
}

particlegrid::~particlegrid(){
	delete grid;
}

void particlegrid::init(const int& x, const int& y, const int& z){
	dimensions = vec3(x,y,z);
	grid = new intgrid(-1);
}

vector<particle*> particlegrid::getCellNeighbors(vec3 index, vec3 numberOfNeighbors){
	//loop through neighbors, for each neighbor, check if cell has particles and push back contents
	vector<particle*> neighbors;
	for( int sx=index.x-numberOfNeighbors.x; sx<=index.x+numberOfNeighbors.x; sx++ ){
		for( int sy=index.y-numberOfNeighbors.y; sy<=index.y+numberOfNeighbors.y; sy++ ) {
			for( int sz=index.z-numberOfNeighbors.z; sz<=index.z+numberOfNeighbors.z; sz++ ) {
				if( sx < 0 || sx > dimensions.x-1 || sy < 0 || sy > dimensions.y-1 || 
					sz < 0 || sz > dimensions.z-1 ){
					continue;
				}
				int cellindex = grid->getCell(vec3(sx, sy, sz));
				if(cellindex>=0){
					int cellparticlecount = cells[cellindex].size();
					for(int a=0; a<cellparticlecount; a++){
						neighbors.push_back(cells[cellindex][a]);
					}
				}	
			}
		}
	}
	// cout << neighbors.size() << endl;
	// cout << index.x << " " << index.y << " " << index.z << endl;
	return neighbors;
}

void particlegrid::sort(const vector<particle*>& particles){
	//clear existing cells
	int cellcount = cells.size();
	for(int i=0; i<cellcount; i++){
		cells[i].clear();
	}

	int particlecount = particles.size();
	int cellscount = cells.size();
	for(int i=0; i<particlecount; i++){
		particle* p = particles[i];
		vec3 pos = p->p;
		pos.x = (int)fmax(0, fmin((int)dimensions.x-1, (int)dimensions.x*pos.x));
		pos.y = (int)fmax(0, fmin((int)dimensions.y-1, (int)dimensions.y*pos.y));
		pos.z = (int)fmax(0, fmin((int)dimensions.z-1, (int)dimensions.z*pos.z));
		int cellindex = grid->getCell(pos);
		// cout << p->p.x << " " << p->p.y << " " << p->p.z << endl;

		if(cellindex>=0){ //if grid has value here, a cell already exists for it
			cells[cellindex].push_back(p);
		}else{ //if grid has no value, create new cell and push index to grid
			vector<particle*> cell;
			cell.push_back(p);
			cells.push_back(cell);
			grid->setCell(pos, cellscount);
			cellscount++;
		}
	}

	// for(int i=0; i<cells.size(); i++){
	// 	cout << cells[i].size() << endl;
	// }
}
