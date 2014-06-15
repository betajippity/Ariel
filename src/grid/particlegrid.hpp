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

class particlegrid{
	public:
		//Initializers
		particlegrid(const glm::vec3& dimensions, const gridtype& type);
		particlegrid(const int& x, const int& y, const int& z, const gridtype& type);
		~particlegrid();

		//Sorting tools
		void sort(std::vector<particle*>& particles);
		std::vector<particle*> getCellNeighbors(glm::vec3 index, glm::vec3 numberOfNeighbors);
		std::vector<particle*> getWallNeighbors(glm::vec3 index, glm::vec3 numberOfNeighbors);

		void markCellTypes(std::vector<particle*>& particles, intgrid* A, float density);
		float cellSDF(int i, int j, int k, float density, geomtype type);

		void buildSDF(macgrid& mgrid, float density);

	private:
		void init(const int& x, const int& y, const int& z, const gridtype& type);

		glm::vec3 dimensions;
		intgrid* grid;
		std::vector< std::vector<particle*> > cells;
		
};
}

#endif
