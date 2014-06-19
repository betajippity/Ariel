// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegrid.cpp
// Helper class for projecting particles back onto the grid and storing particles in cells

#ifndef FLUIDGRID_HPP
#define FLUIDGRID_HPP

#include <tbb/tbb.h>
#include "../utilities/utilities.h"
#include "macgrid.inl"
#include "gridutils.inl"

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class ParticleGrid{
	public:
		//Initializers
		ParticleGrid(const glm::vec3& dimensions);
		ParticleGrid(const int& x, const int& y, const int& z);
		~ParticleGrid();

		//Sorting tools
		void Sort(std::vector<particle*>& particles);
		std::vector<particle*> GetCellNeighbors(const glm::vec3& index, 
												const glm::vec3& numberOfNeighbors);
		std::vector<particle*> GetWallNeighbors(const glm::vec3& index, 
												const glm::vec3& numberOfNeighbors);

		void MarkCellTypes(std::vector<particle*>& particles, Grid<int>* A, const float& density);
		float CellSDF(const int& i, const int& j, const int& k, const float& density, 
					  const geomtype& type);

		void BuildSDF(macgrid& mgrid, const float& density);

	private:
		void Init(const int& x, const int& y, const int& z);

		glm::vec3									m_dimensions;
		Grid<int>*									m_grid;
		std::vector< std::vector<particle*> >		m_cells;
		
};
}

#endif
