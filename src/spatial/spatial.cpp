// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: spatial.cpp. Adapted from Takua Render.
// Implements spatial.hpp

#include "spatial.hpp"

namespace spaceCore {

TraverseAccumulator::TraverseAccumulator(){
    
}

TraverseAccumulator::~TraverseAccumulator(){

}

//The default traverse accumulator does not do anything even remotely interesting
void TraverseAccumulator::RecordIntersection(const rayCore::Intersection& intersect,
                                             const unsigned int& nodeid){
    m_intersection = intersect;
    m_nodeid = nodeid;
}

void TraverseAccumulator::Transform(const glm::mat4& m){
    m_intersection = m_intersection.Transform(m);
}

DebugTraverseAccumulator::DebugTraverseAccumulator(){
    
}

DebugTraverseAccumulator::~DebugTraverseAccumulator(){
    m_intersections.clear();
    m_nodeids.clear();
}

//The default traverse accumulator does not do anything even remotely interesting
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
