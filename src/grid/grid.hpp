// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: grid.hpp
// Templatized grid class

#ifndef GRID_HPP
#define GRID_HPP

#include <tbb/tbb.h>
#include "../utilities/utilities.h"

namespace fluidCore {
//====================================
// Class Declarations
//====================================

template <typename T> class grid{
	public:
		//Initializers
		grid(const glm::vec3& dimensions, const T& background);
		~grid();

		//Cell accessors and setters and whatever
		T getCell(const glm::vec3& index);
		T getCell(const int& x, const int& y, const int& z);

		void setCell(const glm::vec3& index, const T& value);
		void setCell(const int& x, const int& y, const int& z, const T& value);

		void clear();

	protected:
		T*** rawgrid;
		T background;

		glm::vec3 dimensions;
};
}

#include "grid.inl"

#endif
