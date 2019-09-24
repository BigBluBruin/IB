#include "Information_Bottleneck/IB_kernel.h"

IB_kernel::IB_kernel(std::vector<std::vector<double>> input, int quan, int Max_run)
{
    prob_join_xy=input;
    quan_size=quan;
    max_run=Max_run;
}

void IB_kernel::smIB()
{

}