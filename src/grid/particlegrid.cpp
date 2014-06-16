// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegrid.cpp
// Implements particlegrid.hpp

#include "particlegrid.hpp"

namespace fluidCore{

particlegrid::particlegrid(const glm::vec3& dim){
	init((int)dim.x, (int)dim.y, (int)dim.z);
}

particlegrid::particlegrid(const int& x, const int& y, const int& z){
	init(x, y, z);
}

particlegrid::~particlegrid(){
	delete m_grid;
}

void particlegrid::init(const int& x, const int& y, const int& z){
	m_dimensions = glm::vec3(x,y,z);
	m_grid = new grid<int>(glm::vec3(x,y,z), -1);
}

std::vector<particle*> particlegrid::getCellNeighbors(glm::vec3 index,
													  glm::vec3 numberOfNeighbors){
	//loop through neighbors, for each neighbor, check if cell has particles and push back contents
	std::vector<particle*> neighbors;
	for( int sx=index.x-numberOfNeighbors.x; sx<=index.x+numberOfNeighbors.x; sx++ ){
		for( int sy=index.y-numberOfNeighbors.y; sy<=index.y+numberOfNeighbors.y; sy++ ) {
			for( int sz=index.z-numberOfNeighbors.z; sz<=index.z+numberOfNeighbors.z; sz++ ) {
				if( sx < 0 || sx > m_dimensions.x-1 || sy < 0 || sy > m_dimensions.y-1 || 
					sz < 0 || sz > m_dimensions.z-1 ){
					continue;
				}
				int cellindex = m_grid->getCell(glm::vec3(sx, sy, sz));
				if(cellindex>=0){
					int cellparticlecount = m_cells[cellindex].size();
					for(int a=0; a<cellparticlecount; a++){
						neighbors.push_back(m_cells[cellindex][a]);
					}
				}	
			}
		}
	}
	return neighbors;
}

std::vector<particle*> particlegrid::getWallNeighbors(glm::vec3 index, 
													  glm::vec3 numberOfNeighbors){
	std::vector<particle*> neighbors;
	for( int sx=index.x-numberOfNeighbors.x; sx<=index.x+numberOfNeighbors.x-1; sx++ ){
		for( int sy=index.y-numberOfNeighbors.y; sy<=index.y+numberOfNeighbors.y-1; sy++ ) {
			for( int sz=index.z-numberOfNeighbors.z; sz<=index.z+numberOfNeighbors.z-1; sz++ ) {
				if( sx < 0 || sx > m_dimensions.x-1 || sy < 0 || sy > m_dimensions.y-1 || 
					sz < 0 || sz > m_dimensions.z-1 ){
					continue;
				}
				int cellindex = m_grid->getCell(glm::vec3(sx, sy, sz));
				if(cellindex>=0){
					int cellparticlecount = m_cells[cellindex].size();
					for(int a=0; a<cellparticlecount; a++){
						neighbors.push_back(m_cells[cellindex][a]);
					}
				}	
			}
		}
	}
	return neighbors;
}

float particlegrid::cellSDF(int i, int j, int k, float density, geomtype type){
	float accm = 0.0f;
	int cellindex = m_grid->getCell(i,j,k);
	if(cellindex>=0){
		for( int a=0; a<m_cells[cellindex].size(); a++ ) { 
			if( m_cells[cellindex][a]->type == type) {
				accm += m_cells[cellindex][a]->density;
			} else {
				return 1.0f;
			}
		}
	}
	float n0 = 1.0f/(density*density*density);
	return 0.2f*n0-accm;
}

void particlegrid::buildSDF(macgrid& mgrid, float density){
	int x = m_dimensions.x; int y = m_dimensions.y; int z = m_dimensions.z;
	mgrid.L->clear();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				for(int j = 0; j < y; ++j){
					for(int k = 0; k < z; ++k){
						mgrid.L->setCell(i, j, k, cellSDF(i, j, k, density, FLUID));
					}
				}
			}	
		}
	);
}

void particlegrid::markCellTypes(std::vector<particle*>& particles, grid<int>* A, float density){
	int x = m_dimensions.x; int y = m_dimensions.y; int z = m_dimensions.z;
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
				for(int j = 0; j < y; ++j){
					for(int k = 0; k < z; ++k){
						A->setCell(i,j,k, AIR);
						int cellindex = m_grid->getCell(i,j,k);
						if(cellindex>=0 && cellindex<m_cells.size()){
							for( int a=0; a<m_cells[cellindex].size(); a++ ) { 
								if( m_cells[cellindex][a]->type == SOLID ) {
									A->setCell(i,j,k, SOLID);
								}
							}
						}
						if( A->getCell(i,j,k) != SOLID ){
							bool isfluid = cellSDF(i, j, k, density, FLUID) < 0.0 ;
							if(isfluid){
								A->setCell(i,j,k, FLUID);
							}else{
								A->setCell(i,j,k, AIR);
							}
						}
					}
				}
			}
		}
	);
}

void particlegrid::sort(std::vector<particle*>& particles){
	// clear existing cells
	int cellcount = m_cells.size();
	for(int i=0; i<cellcount; i++){
		m_cells[i].clear();
	}

	float maxd = glm::max(glm::max(m_dimensions.x, m_dimensions.y), m_dimensions.z);

	int particlecount = particles.size();
	int cellscount = m_cells.size();
	// cout << particlecount << endl;
	for(int i=0; i<particlecount; i++){
		particle* p = particles[i];

		glm::vec3 pos = p->p;
		pos.x = (int)glm::max(0.0f, glm::min((int)maxd-1.0f, int(maxd)*pos.x));
		pos.y = (int)glm::max(0.0f, glm::min((int)maxd-1.0f, int(maxd)*pos.y));
		pos.z = (int)glm::max(0.0f, glm::min((int)maxd-1.0f, int(maxd)*pos.z));

		int cellindex = m_grid->getCell(pos);
	
		if(cellindex>=0){ //if grid has value here, a cell already exists for it
			m_cells[cellindex].push_back(p);
		}else{ //if grid has no value, create new cell and push index to grid
			std::vector<particle*> cell;
			cell.push_back(p);
			m_cells.push_back(cell);
			m_grid->setCell(pos, cellscount);
			cellscount++;
		}
	}
}
}
