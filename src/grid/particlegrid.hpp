// Kai: FLIP Fluid Simulator
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
		particlegrid(const vec3& dimensions);
		particlegrid(const int& x, const int& y, const int& z);
		~particlegrid();

		//Sorting tools
		void sort(const vector<particle*>& particles);
		vector<particle*> getCellNeighbors(vec3 index, vec3 numberOfNeighbors);

	private:
		void init(const int& x, const int& y, const int& z);

		vec3 dimensions;
		intgrid* grid;
		vector< vector<particle*> > cells;
		
};
}

#endif