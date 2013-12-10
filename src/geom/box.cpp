// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: box.cpp
// Implements box.hpp

#include "box.hpp"

using namespace geomCore;

box::box(const vec3& lowerCorner, const vec3& upperCorner, const geomtype& type){
	this->lowerCorner = lowerCorner;
	this->upperCorner = upperCorner;
	this->type = type;
}

box::box(){
	this->lowerCorner = vec3(0,0,0);
	this->upperCorner = vec3(0,0,0);
	type = SOLID;
}

box::~box(){

}

bool box::isPointInside(const vec3& point, const float& scale){
	if( point.x > (lowerCorner.x*scale) && point.x < (upperCorner.x*scale) &&
	    point.y > (lowerCorner.y*scale) && point.y < (upperCorner.y*scale) &&
	    point.z > (lowerCorner.z*scale) && point.z < (upperCorner.z*scale) ){
		return true;
	}else{
		return false;
	}
}

geomtype box::getType(){
	return type;
}

bool box::isPointInsideWithThickness(const vec3& point, const float& thickness, const float& scale){
	if( point.x > (lowerCorner.x*scale)+thickness && point.x < (upperCorner.x*scale)-thickness &&
	    point.y > (lowerCorner.y*scale)+thickness && point.y < (upperCorner.y*scale)-thickness &&
	    point.z > (lowerCorner.z*scale)+thickness && point.z < (upperCorner.z*scale)-thickness ){
		return true;
	}else{
		return false;
	}
}
