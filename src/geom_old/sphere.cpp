// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: sphere.cpp
// Implements sphere.hpp

#include "sphere.hpp"

using namespace geomCore;

sphere::sphere(const float& radius, const vec3& center, const geomtype& type){
	this->radius = radius;
	this->center = center;
	this->type = type;
}

sphere::sphere(){
	radius = 1.0f;
	center = vec3(0,0,0);
	type = SOLID;
}

sphere::~sphere(){

}

bool sphere::isPointInside(const vec3& point, const float& scale){
	float distanceFromCenter = length(point-center*scale);
	if(distanceFromCenter<radius){
		return true;
	}else{
		return false;
	}
}

geomtype sphere::getType(){
	return type;
}

bool sphere::isPointInsideWithThickness(const vec3& point, const float& thickness, const float& scale){
	float distanceFromCenter = length(point-center*scale);
	if(distanceFromCenter<radius-thickness){
		return true;
	}else{
		return false;
	}	
}
