#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include "Information_Bottleneck/overloadvec.h"
#include "Information_Bottleneck/itbox.h"

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



/*This package gives statistical tools*/
double qfunc(double in);
std::vector<std::vector<double>> gaussian_disretization (double min, double max, int cardi, double sigma2);
std::vector<std::vector<double>> gaussian_disretization2 (double min, double max, int cardi, double sigma21, double sigma22);
std::vector<std::vector<double>> llr_permutation(std::vector<std::vector<double>> &joint_prob, double permutation_factor);
double sum_two_vector(std::vector<std::vector<double>> & input);

