// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: grid.
// Implements grid.hpp

#include "gridutils.inl"

namespace fluidCore{

template <typename T> Grid<T>::Grid(const glm::vec3& dimensions, const T& background){
	m_dimensions = dimensions;
	m_background = background;
	m_rawgrid = CreateGrid<T>(m_dimensions.x+1, m_dimensions.y+1, m_dimensions.z+1);
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)m_dimensions.x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
				for(unsigned int j=0; j<(int)m_dimensions.y+1; ++j){
					for(unsigned int k=0; k<(int)m_dimensions.z+1; ++k){
						m_rawgrid[i][j][k] = m_background;
					}
				}
			}
		}
	);
}

template <typename T> Grid<T>::~Grid(){
	DeleteGrid<T>(m_rawgrid, m_dimensions.x+1, m_dimensions.y+1, m_dimensions.z+1);
}

template <typename T> T Grid<T>::GetCell(const glm::vec3& index){
	return GetCell((int)index.x, (int)index.y, (int)index.z);
}

template <typename T> T Grid<T>::GetCell(const int& x, const int& y, const int& z){
	T cell = m_rawgrid[x][y][z];
	return cell;
}

template <typename T> void Grid<T>::SetCell(const glm::vec3& index, const T& value){
	SetCell((int)index.x, (int)index.y, (int)index.z, value);
}

template <typename T> void Grid<T>::SetCell(const int& x, const int& y, const int& z, 
											const T& value){
	m_rawgrid[x][y][z] = value;
}

template <typename T> void Grid<T>::Clear(){
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)m_dimensions.x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
				for(unsigned int j=0; j<(int)m_dimensions.y+1; ++j){
					for(unsigned int k=0; k<(int)m_dimensions.z+1; ++k){
						m_rawgrid[i][j][k] = m_background;
					}
				}
			}
		}
	);	
}
}

