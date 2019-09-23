#pragma once
#include <algorithm>
#include <functional>
#include <vector>
#include <iterator>
#include <assert.h>



template <typename T>
std::vector<T> operator + (const std::vector<T>& a, const std::vector<T>& b)
{


    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <class T, class Q>
std::vector <T> operator * (const Q c, std::vector <T> &A)
{
    std::vector<T> result;
    result.reserve(A.size());
    for(const auto & a: A)
    {
        result.push_back(c*a);
    }
    return result ;
}

/*template <typename T>
std::vector<T> operator + (const std::vector<T>& a, const std::vector<T>& b);

template <class T, class Q>
std::vector <T> operator * (const Q c, std::vector <T> A);*/



/*template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <class T, class Q>
std::vector <T> operator* (const Q c, std::vector <T> A)
{
    std::transform (A.begin (), A.end (), A.begin (),
                 std::bind1st (std::multiplies <T> () , c)) ;
    return A ;
}*/