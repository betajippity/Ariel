// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: intgrid.hpp
// Class that wraps OpenVDB for underlying grid representation

#ifndef INTGRID_HPP
#define INTGRID_HPP

#include "../utilities/utilities.h"
#include "gridutils.inl"
#include <openvdb/openvdb.h>

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class intgrid{
	public:
		//Initializers
		intgrid(const int& background);
		~intgrid();

		//Cell accessors and setters and whatever
		int getCell(const vec3& index);
		int getCell(const int& x, const int& y, const int& z);

		void setCell(const vec3& index, const int& value);
		void setCell(const int& x, const int& y, const int& z, const int& value);

		openvdb::Int32Grid::Ptr& getVDBGrid();

		void writeVDBGridToFile(string filename);

	private:
		openvdb::Int32Grid::Ptr grid;
		
};
}

#endif
