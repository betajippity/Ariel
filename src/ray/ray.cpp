// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: ray.cpp. Adapted from Takua Render.
// Implements ray.hpp

#include "ray.hpp"

namespace rayCore {

//====================================
// Ray Class
//====================================

HOST DEVICE Ray::Ray(){
    SetContents(glm::vec3(0), glm::vec3(0,1,0), 0.0f, 0);
}

HOST DEVICE Ray::Ray(const glm::vec3& origin, const glm::vec3& direction, const float& frame,
                     const unsigned int& trackingID){
    SetContents(origin, direction, 0.0f, trackingID);
}

HOST DEVICE Ray::Ray(const glm::vec3& origin, const glm::vec3& direction, const float& frame){
    SetContents(origin, direction, 0.0f, 0);
}

HOST DEVICE Ray::~Ray(){

}

HOST DEVICE void Ray::SetContents(const glm::vec3& origin, const glm::vec3& direction, 
                                  const float& frame, const unsigned int& trackingID){
    m_origin = origin;
    m_direction = direction/glm::length(direction);
    m_trackingID = trackingID;
    m_frame = frame;
}

HOST DEVICE glm::vec3 Ray::GetPointAlongRay(const float& distance) const{
    return m_origin + distance * m_direction;
}

HOST DEVICE Ray Ray::Transform(const glm::mat4& m) const{
    glm::vec3 transformedOrigin = glm::vec3(utilityCore::multiply(m, glm::vec4(m_origin, 1.0f)));
    glm::vec3 transformedPoint = glm::vec3(utilityCore::multiply(m, glm::vec4(m_origin+m_direction,
                                                                              1.0f)));
    glm::vec3 transformedDirection = transformedPoint - transformedOrigin;
    return Ray(transformedOrigin, transformedDirection, m_frame, m_trackingID);
}

//====================================
// Intersection Class
//====================================

HOST DEVICE Intersection::Intersection(){
    SetContents(false, glm::vec3(0.0f), glm::vec3(1.0f,0.0f,0.0f), glm::vec2(0.0f), 0, 0);
}

HOST DEVICE Intersection::Intersection(const bool& hit, const glm::vec3& point, 
                                       const glm::vec3& normal, const glm::vec2& uv, 
                                       const unsigned int& objectID, const unsigned int& primID){
    SetContents(hit, point, normal, uv, objectID, primID);
}

HOST DEVICE Intersection::~Intersection(){

}

HOST DEVICE void Intersection::SetContents(const bool& hit, const glm::vec3& point, 
                                           const glm::vec3& normal, const glm::vec2& uv, 
                                           const unsigned int& objectID, 
                                           const unsigned int& primID){
    m_point = point;
    m_hit = hit;
    m_normal = normal;
    m_uv = uv;
    m_objectID = objectID;
    m_primID = primID;
}


HOST DEVICE Intersection Intersection::CompareClosestAgainst(const Intersection& b, 
                                                             const glm::vec3& point){
    if(m_hit && !b.m_hit){
        return *this;
    }else if(!m_hit && b.m_hit){
        return b;
    }else{
        float distanceA = glm::length(m_point - point);
        float distanceB = glm::length(b.m_point - point);
        if(distanceA<=distanceB){
            return *this;
        }else{
            return b;
        }
    }
}

HOST DEVICE Intersection Intersection::Transform(const glm::mat4& m){
    glm::vec3 transformedPoint = glm::vec3(m*glm::vec4(m_point, 1.0f));
    glm::vec3 transformedNormal = glm::vec3(m*glm::vec4(m_normal, 0.0f));
    return Intersection(m_hit, transformedPoint, transformedNormal, m_uv, m_objectID, m_primID); 
}
}
