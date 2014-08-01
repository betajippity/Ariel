// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: geom.cpp
// Implements geom.hpp

#include "geom.hpp"

namespace geomCore {

//====================================
// Geom Class
//====================================

HOST DEVICE Geom::Geom(){
    SetContents(NULL);
}

HOST DEVICE Geom::Geom(GeomInterface* geom){
    SetContents(geom);
}

HOST DEVICE Geom::~Geom(){

}

HOST DEVICE void Geom::SetContents(GeomInterface* geom){
    m_geom = geom;
}

HOST DEVICE void Geom::Intersect(const rayCore::Ray& r, spaceCore::TraverseAccumulator& result){
    m_geom->Intersect(r, result);
}

HOST DEVICE GeomType Geom::GetType(){
    return m_geom->GetType();
}

//====================================
// GeomTransform Class
//====================================

GeomTransform::GeomTransform(){
    SetContents(glm::vec3(0), glm::vec3(0), glm::vec3(1));
}

GeomTransform::GeomTransform(const glm::vec3& t, const glm::vec3& r, const glm::vec3& s){
    SetContents(t, r, s);
}

GeomTransform::~GeomTransform(){

}

void GeomTransform::SetContents(const glm::vec3& t, const glm::vec3& r, const glm::vec3& s){
    m_translation = t;
    m_rotation = r;
    m_scale = s;
}
}
