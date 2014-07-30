// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: shortstack.inl
// Implements shortstack from datastructures.hpp

#ifndef SHORTSTACK_INL
#define SHORTSTACK_INL

#include "datastructures.hpp"

template <typename T> ShortStack<T>::ShortStack(){
    m_currentIndex = -1;
} 

template <typename T> ShortStack<T>::~ShortStack(){

}

template <typename T> void ShortStack<T>::Push(const T& item){
    if(m_currentIndex+1<30){
        m_currentIndex++;
        m_stack[m_currentIndex] = item;
    }
}

template <typename T> T ShortStack<T>::Pop(){
    if(m_currentIndex<0){
        return m_stack[0];
    }else{
        m_currentIndex = m_currentIndex-1;
        return m_stack[m_currentIndex+1];
    }
}

template <typename T> unsigned int ShortStack<T>::Size(){
    return m_currentIndex+1;
}

template <typename T> bool ShortStack<T>::Empty(){
    if(m_currentIndex<0){
        return true;
    }else{
        return false;
    }
}

template <typename T> bool ShortStack<T>::Full(){
    if(m_currentIndex==29){
        return true;
    }else{
        return false;
    }
}

#endif
