// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegrid.cpp
// Helper class for projecting particles back onto the grid and storing what particles are in what cells

#ifndef FLUIDGRID_HPP
#define FLUIDGRID_HPP

#include "../utilities/utilities.h"
#include "macgrid.inl"
#include "gridutils.inl"

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class particlegrid{
	public:
		//Initializers
		particlegrid(const vec3& dimensions, const gridtype& type);
		particlegrid(const int& x, const int& y, const int& z, const gridtype& type);
		~particlegrid();

		//Sorting tools
		void sort(vector<particle*>& particles);
		vector<particle*> getCellNeighbors(vec3 index, vec3 numberOfNeighbors);
		vector<particle*> getWallNeighbors(vec3 index, vec3 numberOfNeighbors);

		void markCellTypes(vector<particle*>& particles, intgrid* A, float density);
		float cellSDF(int i, int j, int k, float density, geomtype type);

		void buildSDF(macgrid& mgrid, float density);

	private:
		void init(const int& x, const int& y, const int& z, const gridtype& type);

		vec3 dimensions;
		intgrid* grid;
		vector< vector<particle*> > cells;
		
};
}

#endif
