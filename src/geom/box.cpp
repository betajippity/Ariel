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

bool box::isPointInside(const vec3& point){
	if( point.x > lowerCorner.x && point.x < upperCorner.x &&
	    point.y > lowerCorner.y && point.y < upperCorner.y &&
	    point.z > lowerCorner.z && point.z < upperCorner.z ){
		return true;
	}else{
		return false;
	}
}

geomtype box::getType(){
	return type;
}

bool box::isPointInsideWithThickness(const vec3& point, const float& thickness){
	if( point.x > lowerCorner.x+thickness && point.x < upperCorner.x-thickness &&
	    point.y > lowerCorner.y+thickness && point.y < upperCorner.y-thickness &&
	    point.z > lowerCorner.z+thickness && point.z < upperCorner.z-thickness ){
		return true;
	}else{
		return false;
	}
}
