// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: perspcam.cpp
// Implements perspcam.hpp

#include "perspcam.hpp"

namespace cameraCore {
    
HOST DEVICE PerspectiveCamera::PerspectiveCamera(){
    SetContents(glm::uvec2(0.0f), glm::vec2(45.0f), 0, 0, 0, NULL);
}

HOST DEVICE PerspectiveCamera::PerspectiveCamera(const glm::uvec2& resolution, 
                                                 const glm::vec2& fov,
                                                 const unsigned int& iterations,
                                                 const unsigned int& traceDepth,
                                                 const unsigned int& numberOfFrames,
                                                 CameraFrame** camFrames){
    SetContents(resolution, fov, iterations, traceDepth, numberOfFrames,
                camFrames);
}

HOST DEVICE PerspectiveCamera::PerspectiveCamera(PerspectiveCameraData data){
    m_id = data.m_id;
    SetContents(data.m_resolution, data.m_fov, data.m_iterations, data.m_traceDepth, 
                data.m_numberOfFrames, data.m_camFrames);
}

HOST DEVICE PerspectiveCamera::~PerspectiveCamera(){

}

HOST DEVICE void PerspectiveCamera::SetContents(const glm::uvec2& resolution, const glm::vec2& fov,
                                                const unsigned int& iterations,
                                                const unsigned int& traceDepth,
                                                const unsigned int& numberOfFrames,
                                                CameraFrame** camFrames){
    m_resolution = resolution;
    m_fov = fov;
    m_iterations = iterations;
    m_traceDepth = traceDepth;
    m_numberOfFrames = numberOfFrames;
    m_camFrames = camFrames;
}
    
HOST DEVICE float PerspectiveCamera::GetAperture(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_aperture;
}

HOST DEVICE float PerspectiveCamera::GetFocal(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_focal;
}

HOST DEVICE glm::vec3 PerspectiveCamera::GetTranslation(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_translation;
}

HOST DEVICE glm::vec3 PerspectiveCamera::GetRotation(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_rotation;
}

HOST DEVICE glm::vec3 PerspectiveCamera::GetView(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_view;
}

HOST DEVICE glm::vec3 PerspectiveCamera::GetUp(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_up;
}

HOST DEVICE float PerspectiveCamera::GetLookat(const unsigned int& frame){
    unsigned int clampedFrame = glm::min(frame, m_numberOfFrames-1);
    return m_camFrames[clampedFrame]->m_lookat;
}

HOST DEVICE CameraType PerspectiveCamera::GetType(){
    return PERSP;
}

HOST DEVICE unsigned int PerspectiveCamera::GetID(){
    return m_id;
}

HOST DEVICE glm::uvec2 PerspectiveCamera::GetResolution(){
    return m_resolution;
}
        
HOST DEVICE glm::vec2 PerspectiveCamera::GetFOV(){
    return m_fov;
}
        
HOST DEVICE unsigned int PerspectiveCamera::GetIterations(){
    return m_iterations;
}

HOST DEVICE unsigned int PerspectiveCamera::GetTraceDepth(){
    return m_traceDepth;
}

HOST DEVICE rayCore::Ray PerspectiveCamera::Raycast(const glm::vec2& xy,
                                                    const glm::vec4& randomSample,
                                                    float& frame){
    int res = glm::min(m_resolution.x, m_resolution.y);
    float fov = glm::min(m_fov.x, m_fov.y);
    int x = glm::floor(res/2.0f)-glm::floor(m_resolution.x/2.0f)+xy.x;
    int y = glm::floor(res/2.0f)-glm::floor(m_resolution.y/2.0f)+xy.y;
    glm::vec3 view = GetView(frame);
    glm::vec3 up = GetUp(frame);
    glm::vec3 P = GetTranslation(frame);
    glm::vec3 horizontalAxis = glm::cross(view, up);
    glm::vec3 verticalAxis = glm::cross(horizontalAxis, view);
    //Get point on our image plane
    glm::vec3 middle = P + view;
    glm::vec3 horizontal = horizontalAxis * glm::tan(fov * 0.5f * float(PI/180));
    glm::vec3 vertical = verticalAxis * glm::tan(-fov * 0.5f * float(PI/180));
    //Convert pixel to normalized device coordinates
    float sx = (x+randomSample.z)/float(res-1);
    float sy = (y+randomSample.w)/float(res-1);
    //Figure out point on image plane based on focal distance
    glm::vec3 pointOnPlaneOneUnitAwayFromEye = middle + (((2.0f*sx)-1.0f) * horizontal) +
    (((2.0f*sy)-1.0f) * vertical);
    glm::vec3 pointOnImagePlane = P + ((pointOnPlaneOneUnitAwayFromEye - P) * GetFocal(frame));
    //Figure out point on lens using random buffer
    glm::vec3 aperturePoint;
    float aperture = GetAperture(frame);
    if(aperture > 0.00001f){ // The small number is an epsilon value.
        //Sample a point on the circular aperture
        float angle = TWO_PI * randomSample.x;
        float distance = aperture * glm::sqrt(randomSample.y);
        float apertureX = glm::cos(angle) * distance;
        float apertureY = glm::sin(angle) * distance;
        aperturePoint = P + (apertureX * horizontalAxis) + (apertureY * verticalAxis);
    }else{
        aperturePoint = P;
    }
    glm::vec3 apertureToImagePlane;
    apertureToImagePlane = pointOnImagePlane - aperturePoint;
    return rayCore::Ray(aperturePoint, apertureToImagePlane, frame);
}
}
