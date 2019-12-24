#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include "Information_Bottleneck/itbox.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include <numeric>

class IB_kernel
{
public:
    std::vector<std::vector<double>> prob_join_xt;
    std::vector<double> prob_t;
    std::vector<unsigned> cluster;
    std::vector<double> threshold;
    double mi;

public:
    std::vector<std::vector<double>> prob_join_xy;
    unsigned quan_size;
    int max_run;

public:
    IB_kernel(std::vector<std::vector<double>> input, unsigned quan, int Max_run);
    void smIB();
    void smIB2(std::vector<double> ref_llr);
    std::vector<unsigned> random_cluster(const unsigned total_num, const unsigned quan_size);
    std::vector<std::vector<double>> quantize_to_xt(std::vector<std::vector<double>> & input, std::vector<unsigned> & cluster);
    unsigned find_threshold(unsigned left_most, unsigned right_most);
    void Progressive_MMI();
};