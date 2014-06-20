// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegrid.cpp
// Implements particlegrid.hpp

#include "particlegrid.hpp"

namespace fluidCore{

ParticleGrid::ParticleGrid(const glm::vec3& dim){
	Init((int)dim.x, (int)dim.y, (int)dim.z);
}

ParticleGrid::ParticleGrid(const int& x, const int& y, const int& z){
	Init(x, y, z);
}

ParticleGrid::~ParticleGrid(){
	delete m_grid;
}

void ParticleGrid::Init(const int& x, const int& y, const int& z){
	m_dimensions = glm::vec3(x,y,z);
	m_grid = new Grid<int>(glm::vec3(x,y,z), -1);
}

std::vector<Particle*> ParticleGrid::GetCellNeighbors(const glm::vec3& index,
													  const glm::vec3& numberOfNeighbors){
	//loop through neighbors, for each neighbor, check if cell has particles and push back contents
	std::vector<Particle*> neighbors;
	for( int sx=index.x-numberOfNeighbors.x; sx<=index.x+numberOfNeighbors.x; sx++ ){
		for( int sy=index.y-numberOfNeighbors.y; sy<=index.y+numberOfNeighbors.y; sy++ ) {
			for( int sz=index.z-numberOfNeighbors.z; sz<=index.z+numberOfNeighbors.z; sz++ ) {
				if( sx < 0 || sx > m_dimensions.x-1 || sy < 0 || sy > m_dimensions.y-1 || 
					sz < 0 || sz > m_dimensions.z-1 ){
					continue;
				}
				int cellindex = m_grid->GetCell(glm::vec3(sx, sy, sz));
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

std::vector<Particle*> ParticleGrid::GetWallNeighbors(const glm::vec3& index, 
													  const glm::vec3& numberOfNeighbors){
	std::vector<Particle*> neighbors;
	for( int sx=index.x-numberOfNeighbors.x; sx<=index.x+numberOfNeighbors.x-1; sx++ ){
		for( int sy=index.y-numberOfNeighbors.y; sy<=index.y+numberOfNeighbors.y-1; sy++ ) {
			for( int sz=index.z-numberOfNeighbors.z; sz<=index.z+numberOfNeighbors.z-1; sz++ ) {
				if( sx < 0 || sx > m_dimensions.x-1 || sy < 0 || sy > m_dimensions.y-1 || 
					sz < 0 || sz > m_dimensions.z-1 ){
					continue;
				}
				int cellindex = m_grid->GetCell(glm::vec3(sx, sy, sz));
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

float ParticleGrid::CellSDF(const int& i, const int& j, const int& k, const float& density, 
							const geomtype& type){
	float accm = 0.0f;
	int cellindex = m_grid->GetCell(i,j,k);
	if(cellindex>=0){
		for( int a=0; a<m_cells[cellindex].size(); a++ ) { 
			if( m_cells[cellindex][a]->m_type == type) {
				accm += m_cells[cellindex][a]->m_density;
			} else {
				return 1.0f;
			}
		}
	}
	float n0 = 1.0f/(density*density*density);
	return 0.2f*n0-accm;
}

void ParticleGrid::BuildSDF(MacGrid& mgrid, const float& density){
	int x = m_dimensions.x; int y = m_dimensions.y; int z = m_dimensions.z;
	mgrid.m_L->Clear();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				for(int j = 0; j < y; ++j){
					for(int k = 0; k < z; ++k){
						mgrid.m_L->SetCell(i, j, k, CellSDF(i, j, k, density, FLUID));
					}
				}
			}	
		}
	);
}

void ParticleGrid::MarkCellTypes(std::vector<Particle*>& particles, Grid<int>* A, 
								 const float& density){
	int x = m_dimensions.x; int y = m_dimensions.y; int z = m_dimensions.z;
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
				for(int j = 0; j < y; ++j){
					for(int k = 0; k < z; ++k){
						A->SetCell(i,j,k, AIR);
						int cellindex = m_grid->GetCell(i,j,k);
						if(cellindex>=0 && cellindex<m_cells.size()){
							for( int a=0; a<m_cells[cellindex].size(); a++ ) { 
								if( m_cells[cellindex][a]->m_type == SOLID ) {
									A->SetCell(i,j,k, SOLID);
								}
							}
						}
						if( A->GetCell(i,j,k) != SOLID ){
							bool isfluid = CellSDF(i, j, k, density, FLUID) < 0.0 ;
							if(isfluid){
								A->SetCell(i,j,k, FLUID);
							}else{
								A->SetCell(i,j,k, AIR);
							}
						}
					}
				}
			}
		}
	);
}

void ParticleGrid::Sort(std::vector<Particle*>& particles){
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
		Particle* p = particles[i];

		glm::vec3 pos = p->m_p;
		pos.x = (int)glm::max(0.0f, glm::min((int)maxd-1.0f, int(maxd)*pos.x));
		pos.y = (int)glm::max(0.0f, glm::min((int)maxd-1.0f, int(maxd)*pos.y));
		pos.z = (int)glm::max(0.0f, glm::min((int)maxd-1.0f, int(maxd)*pos.z));

		int cellindex = m_grid->GetCell(pos);
	
		if(cellindex>=0){ //if grid has value here, a cell already exists for it
			m_cells[cellindex].push_back(p);
		}else{ //if grid has no value, create new cell and push index to grid
			std::vector<Particle*> cell;
			cell.push_back(p);
			m_cells.push_back(cell);
			m_grid->SetCell(pos, cellscount);
			cellscount++;
		}
	}
}
}
