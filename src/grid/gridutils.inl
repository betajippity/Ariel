// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: gridutils.inl
// Contains some handy macros for looping over grid cells

#ifndef GRIDUTILS_INL
#define GRIDUTILS_INL

#define FOR_EACH_CELL(x, y, z) \
	for(int i = 0; i < x; i++) \
		for(int j = 0; j < y; j++) \
			for(int k = 0; k < z; k++)

#define FOR_EACH_FLOW_X(x, y, z) \
	for(int i = 0; i < x+1; i++)  \
	  	for(int j = 0; j < y; j++) \
	    	for(int k = 0; k < z; k++) 

#define FOR_EACH_FLOW_Y(x, y, z) \
 	for(int i = 0; i < x; i++)  \
	  	for(int j = 0; j < y+1; j++) \
	    	for(int k = 0; k < z; k++) 

#define FOR_EACH_FLOW_Z(x, y, z) \
	for(int i = 0; i < x; i++)  \
	  	for(int j = 0; j < y; j++) \
	    	for(int k = 0; k < z+1; k++) 

#endif
