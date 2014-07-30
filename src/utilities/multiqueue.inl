// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: multiqueue.inl
// Implements multiqueue from datastructures.hpp

#ifndef MULTIQUEUE_INL
#define MULTIQUEUE_INL

#include "datastructures.hpp"

template <typename T> MultiQueue<T>::MultiQueue(){
    
} 

template <typename T> MultiQueue<T>::~MultiQueue(){
    m_queue.clear();
}

template <typename T> void MultiQueue<T>::Push(const T& item){
    m_queue.push(item);
}

template <typename T> void MultiQueue<T>::Push(const std::vector<T>& items){
    unsigned int itemsCount = items.size();
    for(int i=0; i<itemsCount; i++){
        m_queue.push(items[i]);
    }
}

template <typename T> T MultiQueue<T>::Pop(){
    T item;
    m_queue.try_pop(item);
    return item;
}

template <typename T> std::vector<T> MultiQueue<T>::Pop(const unsigned int& count){
    std::vector<T> items(count);
    for(int i=0; i<count; i++){
        m_queue.try_pop(items[i]);
    }
    return items;
} 

template <typename T> unsigned int MultiQueue<T>::Size(){
    unsigned int size = m_queue.unsafe_size();
    return size; 
}

template <typename T> bool MultiQueue<T>::Empty(){
    bool empty = m_queue.empty();
    return empty;
}

#endif
