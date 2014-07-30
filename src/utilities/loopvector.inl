// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: loopvector.inl
// Implements loopvector from datastructures.hpp

#ifndef LOOPVECTOR_INL
#define LOOPVECTOR_INL

#include "datastructures.hpp"

template <typename T> LoopVector<T>::LoopVector(){
	m_currentIndex = 0;
	m_size = 0;
}

template <typename T> LoopVector<T>::~LoopVector(){
	m_vector.clear();
}

template <typename T> void LoopVector<T>::PushBack(const T& item){
	m_vector.push_back(item);
	m_size = m_size + 1;
}

template <typename T> T LoopVector<T>::GetElement(){
	m_currentIndex.compare_and_swap(0, m_size);
	unsigned int i = m_currentIndex.fetch_and_add(1);
	return m_vector[i];
}

#endif
