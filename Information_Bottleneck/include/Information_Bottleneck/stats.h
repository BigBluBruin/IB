#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include "Information_Bottleneck/overloadvec.h"
#include "Information_Bottleneck/itbox.h"


/*This package gives statistical tools*/
double qfunc(double in);

std::vector<std::vector<double>> gaussian_disretization (double min, double max, int cardi, double sigma2);

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}