// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: spatial.cpp. Adapted from Takua Render.
// Implements spatial.hpp

#include "spatial.hpp"

namespace spaceCore {

TraverseAccumulator::TraverseAccumulator(){
    m_origin = glm::vec3(0.0f);   
}

TraverseAccumulator::TraverseAccumulator(const glm::vec3& origin){
    m_origin = origin;  
}

TraverseAccumulator::~TraverseAccumulator(){

}

//The default traverse accumulator does not do anything even remotely interesting
void TraverseAccumulator::RecordIntersection(const rayCore::Intersection& intersect,
                                             const unsigned int& nodeid){
    if(m_intersection.m_hit==false && intersect.m_hit==true){
        m_intersection = intersect;
        m_nodeid = nodeid;   
    }else if(intersect.m_hit==true){
        float currentDistance = glm::length(m_intersection.m_point - m_origin);
        float newDistance = glm::length(intersect.m_point - m_origin);
        if(newDistance<currentDistance){
            m_intersection = intersect;
            m_nodeid = nodeid;   
        }
    }
}

void TraverseAccumulator::Transform(const glm::mat4& m){
    m_intersection = m_intersection.Transform(m);
    m_origin = glm::vec3(utilityCore::multiply(m, glm::vec4(m_origin, 1.0f)));
}

HitCountTraverseAccumulator::HitCountTraverseAccumulator(){
    m_numberOfHits = 0;
    m_origin = glm::vec3(0.0f);
}

HitCountTraverseAccumulator::HitCountTraverseAccumulator(const glm::vec3& origin){
    m_numberOfHits = 0;
    m_origin = origin;  
}

HitCountTraverseAccumulator::~HitCountTraverseAccumulator(){

}

//Just like the default traverse accumulator, but also tracks the number of
//times a hit has been accumulated. Useful for checking if things have
//passed through a mesh or not
void HitCountTraverseAccumulator::RecordIntersection(const rayCore::Intersection& intersect,
                                                     const unsigned int& nodeid){
    if(m_intersection.m_hit==false && intersect.m_hit==true){
        m_intersection = intersect;
        m_nodeid = nodeid;   
        m_numberOfHits++;
    }else if(intersect.m_hit==true){
        float currentDistance = glm::length(m_intersection.m_point - m_origin);
        float newDistance = glm::length(intersect.m_point - m_origin);
        if(newDistance<currentDistance){
            m_intersection = intersect;
            m_nodeid = nodeid;  
        }
        m_numberOfHits++; 
    }
}

void HitCountTraverseAccumulator::Transform(const glm::mat4& m){
    m_intersection = m_intersection.Transform(m);
    m_origin = glm::vec3(utilityCore::multiply(m, glm::vec4(m_origin, 1.0f)));
}

DebugTraverseAccumulator::DebugTraverseAccumulator(){
    
}

DebugTraverseAccumulator::~DebugTraverseAccumulator(){
    m_intersections.clear();
    m_nodeids.clear();
}

//Keeps a copy of every single hit that is recorded to this accumulator
void DebugTraverseAccumulator::RecordIntersection(const rayCore::Intersection& intersect,
                                                  const unsigned int& nodeid){
    m_intersections.push_back(intersect);
    m_nodeids.push_back(nodeid);

}

void DebugTraverseAccumulator::Transform(const glm::mat4& m){
    unsigned int intersectionCount = m_intersections.size();
    for(unsigned int i=0; i<intersectionCount; i++){
        m_intersections[i] = m_intersections[i].Transform(m);
    }
}
}
