// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: perspcam.hpp
// Perspective camera implements CameraInterface

#ifndef PERSPCAM_HPP
#define PERSPCAM_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include "../camera.hpp"

namespace cameraCore {

//====================================
// Struct Declarations
//====================================

//Only used as a temp data transfer container for CUDA memory setup
struct PerspectiveCameraData {
    glm::uvec2          m_resolution;
    glm::vec2           m_fov;
    unsigned int        m_iterations;
    unsigned int        m_traceDepth;
    unsigned int        m_numberOfFrames;
    CameraFrame**       m_camFrames;
    unsigned int        m_id;
};

//====================================
// Class Declarations
//====================================

class PerspectiveCamera: public virtual CameraInterface {
    public:
        HOST DEVICE PerspectiveCamera();
    HOST DEVICE PerspectiveCamera(const glm::uvec2& resolution, const glm::vec2& fov, 
                                      const unsigned int& iterations, 
                                      const unsigned int& traceDepth, 
                                      const unsigned int& numberOfFrames, CameraFrame** camFrames);
        HOST DEVICE PerspectiveCamera(PerspectiveCameraData data);
        HOST DEVICE ~PerspectiveCamera();

        HOST DEVICE void SetContents(const glm::uvec2& resolution, const glm::vec2& fov,
                                     const unsigned int& iterations, 
                                     const unsigned int& traceDepth,
                                     const unsigned int& numberOfFrames, CameraFrame** camFrames);
        HOST DEVICE rayCore::Ray Raycast(const glm::vec2& xy, const glm::vec4& sample,
                                         float& frame);
        HOST DEVICE float GetAperture(const unsigned int& frame);
        HOST DEVICE float GetFocal(const unsigned int& frame);
        HOST DEVICE glm::vec3 GetTranslation(const unsigned int& frame);
    HOST DEVICE glm::vec3 GetRotation(const unsigned int& frame);
    HOST DEVICE glm::vec3 GetView(const unsigned int& frame);
    HOST DEVICE glm::vec3 GetUp(const unsigned int& frame);
    HOST DEVICE float GetLookat(const unsigned int& frame);
        HOST DEVICE unsigned int GetID();
        HOST DEVICE CameraType GetType();
        HOST DEVICE glm::uvec2 GetResolution();
        HOST DEVICE glm::vec2 GetFOV();
        HOST DEVICE unsigned int GetIterations();
        HOST DEVICE unsigned int GetTraceDepth();
        
        glm::uvec2          m_resolution;
    glm::vec2           m_fov;
    unsigned int        m_iterations;
    unsigned int        m_traceDepth;
    unsigned int        m_numberOfFrames;
    CameraFrame**       m_camFrames;
    unsigned int        m_id;
};
}

#endif
