// UtilityCore: A utility library. Part of the TAKUA Render project.
// Written by Yining Karl Li
// Version 0.5.13.39a
//  
// File: utilities.inl
// A collection/kitchen sink of generally useful functions

#ifndef UTILITIES_INL
#define UTILITIES_INL

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

#include <glm/gtc/matrix_transform.inl>
#include <glm/glm.hpp>
#include <iostream>
#include <sys/timeb.h>
#include <cstdio>
#include <cstring>
#include <fstream>

//====================================
// Math stuff
//====================================

float utilityCore::clamp(float f, float min, float max){
    if(f<min){
        return min;
    }else if(f>max){
        return max;
    }else{
        return f;
    }
}

glm::vec3 utilityCore::clampRGB(glm::vec3 color){
    if(color[0]<0){
        color[0]=0;
    }else if(color[0]>255){
        color[0]=255;
    }
    if(color[1]<0){
        color[1]=0;
    }else if(color[1]>255){
        color[1]=255;
    }
    if(color[2]<0){
        color[2]=0;
    }else if(color[2]>255){
        color[2]=255;
    }
    return color;
}

HOST DEVICE bool utilityCore::epsilonCheck(float a, float b){
    if(glm::abs(glm::abs(a)-glm::abs(b))<EPSILON){
        return true;
    }else{
        return false;
    }
}

//====================================
// String wrangling stuff
//====================================

std::string utilityCore::padString(int length, std::string str){
    int strlength = str.length();
    std::string pad = "";
    for(int i=strlength; i<length; i++){
        pad = pad + "0";
    }
    return pad + str;
}

bool utilityCore::replaceString(std::string& str, const std::string& from, const std::string& to){
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::string utilityCore::convertIntToString(int number){
    std::stringstream ss;
    ss << number;
    return ss.str();
}

std::vector<std::string> utilityCore::tokenizeString(std::string str, std::string separator){
    std::vector<std::string> results;
    char * cstr, *p;
    std::string strt = str;
    cstr = new char[strt.size()+1];
    std::strcpy (cstr, strt.c_str());
    p=std::strtok (cstr, separator.c_str());
    while (p!=NULL){
        results.push_back(p);
        p=strtok(NULL, separator.c_str());
    }
    delete [] cstr;
    delete [] p;
    return results;
}

std::vector<std::string> utilityCore::tokenizeStringByAllWhitespace(std::string str){
    std::stringstream strstr(str);
    std::istream_iterator<std::string> it(strstr);
    std::istream_iterator<std::string> end;
    std::vector<std::string> results(it, end);
    return results;
}

std::string utilityCore::getLastNCharactersOfString(std::string s, int n){
    return s.substr(s.length()-n, n);
}

std::string utilityCore::getFirstNCharactersOfString(std::string s, int n){
    return s.substr(0, n);
}

//====================================
// Time stuff
//====================================

int utilityCore::getMilliseconds(){
    timeb tb;
    ftime( &tb );
    int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
    return nCount;
}

int utilityCore::compareMilliseconds(int referenceTime){
    int elapsedMilliseconds = getMilliseconds() - referenceTime;
    if ( elapsedMilliseconds < 0 ){
        elapsedMilliseconds += 0x100000 * 1000;
    }
    return elapsedMilliseconds;
}

//====================================
// Matrix stuff
//====================================

glm::mat4 utilityCore::buildTransformationMatrix(glm::vec3 translation, glm::vec3 rotation,
                                                 glm::vec3 scale){
    glm::mat4 translationMat = glm::translate(glm::mat4(), translation);
    glm::mat4 rotationMat = glm::rotate(glm::mat4(), rotation.z, glm::vec3(0,0,1));
    rotationMat = rotationMat*glm::rotate(glm::mat4(), rotation.y, glm::vec3(0,1,0));
    rotationMat = rotationMat*glm::rotate(glm::mat4(), rotation.x, glm::vec3(1,0,0));
    glm::mat4 scaleMat = glm::scale(glm::mat4(), scale);
    return translationMat*rotationMat*scaleMat;
}

glm::mat4 utilityCore::buildInverseTransformationMatrix(glm::vec3 translation, glm::vec3 rotation, 
                                                        glm::vec3 scale){
    return glm::inverse(buildTransformationMatrix(translation, rotation, scale));
}

HOST DEVICE glm::vec4 utilityCore::multiply(glm::mat4 m, glm::vec4 v){
    glm::vec4 r;
    r.x = m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z + m[3][0] * v.w;
    r.y = m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z + m[3][1] * v.w;
    r.z = m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z + m[3][2] * v.w;
    r.w = m[0][3] * v.x + m[1][3] * v.y + m[2][3] * v.z + m[3][3] * v.w;
    return r;
}

//====================================
// GLM Printers
//====================================

void utilityCore::printMat4(glm::mat4 m){
    std::cout << m[0][0] << " " << m[1][0] << " " << m[2][0] << " " << m[3][0] << std::endl;
    std::cout << m[0][1] << " " << m[1][1] << " " << m[2][1] << " " << m[3][1] << std::endl;
    std::cout << m[0][2] << " " << m[1][2] << " " << m[2][2] << " " << m[3][2] << std::endl;
    std::cout << m[0][3] << " " << m[1][3] << " " << m[2][3] << " " << m[3][3] << std::endl;
}

void utilityCore::printVec4(glm::vec4 m){
    std::cout << m[0] << " " << m[1] << " " << m[2] << " " << m[3] << std::endl;
}

void utilityCore::printVec3(glm::vec3 m){
    std::cout << m[0] << " " << m[1] << " " << m[2] << std::endl;
}

//====================================
// GL Stuff
//====================================

void utilityCore::fovToPerspective(float fovy, float aspect, float zNear, glm::vec2& xBounds, 
                                   glm::vec2& yBounds){
    yBounds[1] = zNear * tan(fovy*(float)PI/360.0f);
    yBounds[0] = -yBounds[1];
    xBounds[0] = yBounds[0]*aspect;
    xBounds[1] = yBounds[1]*aspect;
}

//====================================
// IO Stuff
//====================================

std::string utilityCore::readFileAsString(std::string filename){
    std::ifstream t(filename.c_str());
    std::string str;
    t.seekg(0, std::ios::end);   
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);
    str.assign((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
    return str;
}

std::string utilityCore::getRelativePath(std::string path){
    std::string relativePath;
    std::vector<std::string> pathTokens = utilityCore::tokenizeString(path, "/");
    for(int i=0; i<pathTokens.size()-1; i++){
        relativePath = relativePath + pathTokens[i] + "/";
    }
    return relativePath;
}

#endif
