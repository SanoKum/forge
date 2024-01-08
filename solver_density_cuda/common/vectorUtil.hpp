#ifndef VECTOR_UTIL_H
#define VECTOR_UTIL_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

template<typename T>
bool ifEqualComponent(const std::vector<T> &v1 , const std::vector<T> &v2 ) 
{
    if (v1.size() != v2.size()) return false;

    for (auto it1 = v1.begin(), e=v1.end() ; it1 != e ; ++it1)
    {
        bool isContained = false;
        for (auto it2 = v2.begin(), e2=v2.end() ; it2 != e2 ; ++it2)
        {
            //cout << *it1 << " " << *it2 << endl;
            if (*it1 == *it2) 
            {
                isContained = true;
                break;
            };
        }
        if (isContained == false) return false;
    }

    return true;
}

template<typename T>
void print1Dvector(const std::vector<T> &v1 ) 
{
    for (auto it = v1.begin(), e=v1.end(); it != e; ++it)
    {
        cout << *it <<  " ";
    }
    cout << endl;
}

template<typename T>
void print2Dvector(const std::vector<std::vector<T>> &v1 ) 
{

    for (auto it = v1.begin(), e=v1.end(); it != e; ++it)
    {
        for (auto it2 = (*it).begin(), e2=(*it).end(); it2 != e2; ++it2)
        {
            cout << *it2 <<  " ";
        }
        cout << endl;
    }
}


template<typename T>
int findIndex(const std::vector<T> &v1 , const T &ele ) 
{
    auto itr = std::find(v1.begin(), v1.end(), ele);

    if (itr == v1.end()) 
    {
        return -1;
    } else {
        int wantedIndex = std::distance(v1.begin(), itr);
        return wantedIndex;
    }
}

#endif