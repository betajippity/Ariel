// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: aabb.cpp. Adapted from Takua Render.
// Implements aabb.hpp

#include "aabb.hpp"

namespace spaceCore {

Aabb::Aabb(){
    SetContents(glm::vec3(REALLY_BIG_NUMBER), glm::vec3(-REALLY_BIG_NUMBER), glm::vec3(0), -1);
}

Aabb::Aabb(const glm::vec3& min, const glm::vec3& max, const int& id){
    glm::vec3 center  = (min + max)/2.0f;
    SetContents(min, max, center, id);
}


Aabb::Aabb(const glm::vec3& min, const glm::vec3& max, const glm::vec3& centroid, const int& id){
    SetContents(min, max, centroid, id);
}

Aabb::~Aabb(){

}

HOST DEVICE float Aabb::FastIntersectionTest(const rayCore::Ray& r){
    if(r.m_origin.x >= m_min.x && r.m_origin.y >= m_min.y && r.m_origin.z >= m_min.z &&
       r.m_origin.x <= m_max.x && r.m_origin.y <= m_max.y && r.m_origin.z <= m_max.z){
        return 0;
    }
    glm::vec3 rdirection = r.m_direction;
    glm::vec3 dirfrac;
    dirfrac.x = 1.0f / rdirection.x;
    dirfrac.y = 1.0f / rdirection.y;
    dirfrac.z = 1.0f / rdirection.z;
    float t1 = (m_min.x - r.m_origin.x)*dirfrac.x;
    float t2 = (m_max.x - r.m_origin.x)*dirfrac.x;
    float t3 = (m_min.y - r.m_origin.y)*dirfrac.y;
    float t4 = (m_max.y - r.m_origin.y)*dirfrac.y;
    float t5 = (m_min.z - r.m_origin.z)*dirfrac.z;
    float t6 = (m_max.z - r.m_origin.z)*dirfrac.z;
    float tmin = glm::max(glm::max(glm::min(t1, t2), glm::min(t3, t4)), glm::min(t5, t6));
    float tmax = glm::min(glm::min(glm::max(t1, t2), glm::max(t3, t4)), glm::max(t5, t6));
    if(tmax < 0){
        return -1;
    }
    if(tmin > tmax){
        return -1;
    }
    return tmin;
}

Aabb Aabb::Transform(const glm::mat4& transform){
    //transform all corners of aabb
    glm::vec3 tmin = glm::vec3(utilityCore::multiply(transform, glm::vec4(m_min, 1.0f)));            
    glm::vec3 tmax = glm::vec3(utilityCore::multiply(transform, glm::vec4(m_max, 1.0f)));            
    glm::vec3 tcentroid = glm::vec3(utilityCore::multiply(transform, 
                                                          glm::vec4(m_centroid, 1.0f)));            
    glm::vec3 m0 = glm::vec3(utilityCore::multiply(transform, 
                                                   glm::vec4(m_max.x, m_min.y, m_min.z, 1.0f)));    
    glm::vec3 m1 = glm::vec3(utilityCore::multiply(transform, 
                                                   glm::vec4(m_max.x, m_min.y, m_max.z, 1.0f)));    
    glm::vec3 m2 = glm::vec3(utilityCore::multiply(transform, 
                                                   glm::vec4(m_max.x, m_max.y, m_min.z, 1.0f)));    
    glm::vec3 m3 = glm::vec3(utilityCore::multiply(transform, 
                                                   glm::vec4(m_min.x, m_min.y, m_max.z, 1.0f)));    
    glm::vec3 m4 = glm::vec3(utilityCore::multiply(transform, 
                                                   glm::vec4(m_min.x, m_max.y, m_min.z, 1.0f)));    
    glm::vec3 m5 = glm::vec3(utilityCore::multiply(transform, 
                                                   glm::vec4(m_min.x, m_max.y, m_max.z, 1.0f)));
    //build new aabb, return
    Aabb a;
    a.m_min = glm::min(glm::min(glm::min(tmin, tmax), glm::min(m0, m1)), 
                       glm::min(glm::min(m2, m3), glm::min(m4, m5)));
    a.m_max = glm::max(glm::max(glm::max(tmin, tmax), glm::max(m0, m1)), 
                       glm::max(glm::max(m2, m3), glm::max(m4, m5)));
    a.m_id = m_id;
    a.m_centroid = tcentroid; 
    return a; 
}

double Aabb::CalculateSurfaceArea(){
	double xlength = m_max.x-m_min.x;
	double ylength = m_max.y-m_min.y;
	double zlength = m_max.z-m_min.z;
	return 2.0f*((xlength*ylength)+(ylength*zlength)+(zlength*xlength));
}

void Aabb::SetContents(const glm::vec3& min, const glm::vec3& max, const glm::vec3& centroid,
                       const int& id){
    m_min = min;
    m_max = max;
    m_centroid = centroid;
    m_id = id;
}

void Aabb::ExpandAabb(const glm::vec3& exMin, const glm::vec3& exMax){
    m_min = glm::min(m_min, exMin);
    m_max = glm::max(m_max, exMax);
}
}

