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
#include <openvdb/tools/Interpolation.h>

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class floatgrid{
	public:
		//Initializers
		floatgrid(const gridtype& type, const vec3& dimensions, const int& background);
		floatgrid();
		floatgrid(openvdb::FloatGrid::Ptr grid);
		~floatgrid();

		//Cell accessors and setters and whatever
		float getCell(const vec3& index);
		float getCell(const int& x, const int& y, const int& z);

		void setCell(const vec3& index, const float& value);
		void setCell(const int& x, const int& y, const int& z, const float& value);

		float getInterpolatedCell(const vec3& index);
		float getInterpolatedCell(const float& x, const float& y, const float& z);

		openvdb::FloatGrid::Ptr& getVDBGrid();

		void writeVDBGridToFile(string filename);

		void clear();

	protected:
		openvdb::FloatGrid::Ptr vdbgrid;
		float*** rawgrid;

		gridtype type;
		vec3 dimensions;

		int background;
};
}

#endif
