// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: intgrid.hpp
// Class that wraps OpenVDB for underlying grid representation

#ifndef INTGRID_HPP
#define INTGRID_HPP

#include <openvdb/openvdb.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include "../utilities/utilities.h"
#include "gridutils.inl"

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class intgrid{
	public:
		//Initializers
		intgrid(const gridtype& type, const glm::vec3& dimensions, const int& background);
		~intgrid();

		//Cell accessors and setters and whatever
		int getCell(const glm::vec3& index);
		int getCell(const int& x, const int& y, const int& z);

		void setCell(const glm::vec3& index, const int& value);
		void setCell(const int& x, const int& y, const int& z, const int& value);

		gridtype getGridType();

		openvdb::Int32Grid::Ptr& getVDBGrid();
		void writeVDBGridToFile(std::string filename);

	private:
		openvdb::Int32Grid::Ptr vdbgrid;
		int*** rawgrid;

		gridtype type;
		glm::vec3 dimensions;

		int background;

		tbb::mutex SetCellLock;
};
}

#endif
