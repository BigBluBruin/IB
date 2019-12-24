#pragma once
#include <algorithm>
#include <functional>
#include <vector>
#include <iterator>
#include <assert.h>
#include <iostream>



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

template <typename T>
std::vector<std::vector<T>> operator + (const std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b)
{
    //verify correction
        std::size_t a_out=a.size();
        std::size_t a_in=a[0].size();
        std::size_t b_out=b.size();
        std::size_t b_in=b[0].size();
        if (a_out!=b_out||a_in!=b_in)
            std::cout<<"joint prob. wrong size: sizes are not compitable ..."<<std::endl;
    //add operation
    std::vector<std::vector<T>> result(a_out, std::vector<T> (a_in,0));
    for (std::size_t out=0; out<a_out; out++)
        for(std::size_t in=0; in<a_in; in++)
            {
                result[out][in]=a[out][in]+b[out][in];
            }
    return result;
}

template <typename T>
double	vectors_distance(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<double>	auxiliary;

	std::transform (a.begin(), a.end(), b.begin(), std::back_inserter(auxiliary),//
	[](T element1, T element2) {return pow((element1-element2),2);});
	auxiliary.shrink_to_fit();

	return  sqrt(std::accumulate(auxiliary.begin(), auxiliary.end(), 0));
} // end template vectors_distance



/*
Consider the following (n,k,d) block code:
D0  D1  D2  D3  D4   | P0
D5  D6  D7  D8  D9   | P1 
D10 D11 D12 D13 D14  | P2 
-------------------------
P3  P4  P5  P6  P7   |
where D0-D14 are data bits, P0-P2 are row parity bits and P3-P7 are column parity bits. The transmitted code word will be:
D0 D1 D2 ... D13 D14 P0 P1 ... P6 P7

How many errors and erasures can it correct? Write a syndrome decoding program that takes a 23 bit input, find the syndrome, determine the error vector from the syndrome and estimate the most likely original codeword.
*/
