/*
This package is used to realize parobability combination
the following parts are includes
    1. probability combination
    2. LLR infusion
    3. sort joint probabilty based on their LLR value
    4. promise symmetry
*/
#pragma once
#include "vector"
#include "iostream"
#include "Information_Bottleneck/itbox.h"
#include <iterator>
#include <algorithm>
#include <numeric>


std::vector<std::vector<double>> prob_combination (std::vector<std::vector<double>> & first_input, std::vector<std::vector<double>> & second_input, const char type[]);
void prob_sort(std::vector<std::vector<double>> & input);
std::vector<unsigned> sorted_pos(std::vector<double> &input);
std::vector<std::vector<double>> prob_approx(std::vector<std::vector<double>> input);
std::vector<std::vector<double>> llr_combination(std::vector<std::vector<double>> & input, double threshold);
std::vector<std::vector<double>> clip_prob(std::vector<std::vector<double>> & input, double threshold);