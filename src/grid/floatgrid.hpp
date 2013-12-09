// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: floatgrid.hpp
// Class that wraps OpenVDB for underlying grid representation

#ifndef FLOATGRID_HPP
#define FLOATGRID_HPP

#include "../utilities/utilities.h"
#include "gridutils.inl"
#include <openvdb/openvdb.h>

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class floatgrid{
	public:
		//Initializers
		floatgrid(const float& background);
		~floatgrid();

		//Cell accessors and setters and whatever
		float getCell(const vec3& index);
		float getCell(const int& x, const int& y, const int& z);

		void setCell(const vec3& index, const float& value);
		void setCell(const int& x, const int& y, const int& z, const float& value);

	private:
		openvdb::FloatGrid::Ptr grid;

};
}

#endif
