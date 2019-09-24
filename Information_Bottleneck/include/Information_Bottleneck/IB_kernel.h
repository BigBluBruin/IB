#pragma once
#include <iostream>
#include <fstream>
#include <vector>

class IB_kernel
{
private:
    std::vector<std::vector<double>> prob_join_xt;
    std::vector<double> prob_t;
    double mi;

public:
    std::vector<std::vector<double>> prob_join_xy;
    int quan_size;
    int max_run;

public:
    IB_kernel(std::vector<std::vector<double>> input, int quan, int Max_run);
    void smIB();
    std::vector<int> random_cluster();
};