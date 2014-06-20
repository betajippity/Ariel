// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: kernels.inl
// Some simple filtering kernels

#ifndef KERNELS_INL
#define KERNELS_INL

#include "../utilities/utilities.h"

namespace mathCore {
//====================================
// Struct and Function Declarations
//====================================

//Forward declarations for externed inlineable methods
extern inline float Smooth(const float& r2, const float& h);
extern inline float Sharpen(const float& r2, const float& h);
extern inline float Sqrlength(const glm::vec3& p0, const glm::vec3& p1);

//====================================
// Function Implementations
//====================================
float Sqrlength(const glm::vec3& p0, const glm::vec3& p1){
	float a = p0.x - p1.x;
	float b = p0.y - p1.y;
	float c = p0.z - p1.z;
	return a*a + b*b + c*c;
}

float Smooth(const float& r2, const float& h) {
    return glm::max(1.0f-r2/(h*h), 0.0f);
}

float Sharpen(const float& r2, const float& h) {
    return glm::max(h*h/glm::max(r2,(float)1.0e-5) - 1.0f, 0.0f);
}
}

#endif
