// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: grid.
// Implements grid.hpp

#include "gridutils.inl"

namespace fluidCore{

template <typename T> grid<T>::grid(const glm::vec3& dimensions, const T& background){
	this->dimensions = dimensions;
	this->background = background;
	rawgrid = createGrid<T>(dimensions.x+1, dimensions.y+1, dimensions.z+1);
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)dimensions.x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
				for(unsigned int j=0; j<(int)dimensions.y+1; ++j){
					for(unsigned int k=0; k<(int)dimensions.z+1; ++k){
						rawgrid[i][j][k] = background;
					}
				}
			}
		}
	);
}

template <typename T> grid<T>::~grid(){
	deleteGrid<T>(rawgrid, dimensions.x+1, dimensions.y+1, dimensions.z+1);
}

template <typename T> T grid<T>::getCell(const glm::vec3& index){
	return getCell((int)index.x, (int)index.y, (int)index.z);
}

template <typename T> T grid<T>::getCell(const int& x, const int& y, const int& z){
	T cell = rawgrid[x][y][z];
	return cell;
}

template <typename T> void grid<T>::setCell(const glm::vec3& index, const T& value){
	setCell((int)index.x, (int)index.y, (int)index.z, value);
}

template <typename T> void grid<T>::setCell(const int& x, const int& y, const int& z, 
											const T& value){
	rawgrid[x][y][z] = value;
}

template <typename T> void grid<T>::clear(){
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(int)dimensions.x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
				for(unsigned int j=0; j<(int)dimensions.y+1; ++j){
					for(unsigned int k=0; k<(int)dimensions.z+1; ++k){
						rawgrid[i][j][k] = background;
					}
				}
			}
		}
	);	
}
}

