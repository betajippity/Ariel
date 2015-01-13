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

template <typename T> class Grid{
    public:
        //Initializers
        Grid(const glm::vec3& dimensions, const T& background);
        ~Grid();

        //Cell accessors and setters and whatever
        T GetCell(const glm::vec3& index);
        T GetCell(const int& x, const int& y, const int& z);

        void SetCell(const glm::vec3& index, const T& value);
        void SetCell(const int& x, const int& y, const int& z, const T& value);

        void Clear();

    protected:
        T***        m_rawgrid;
        T           m_background;

        glm::vec3   m_dimensions;
};
}

#include "grid.inl"

#endif
