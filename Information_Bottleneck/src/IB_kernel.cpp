#include "Information_Bottleneck/IB_kernel.h"

IB_kernel::IB_kernel(std::vector<std::vector<double>> input, unsigned quan, int Max_run)
{
    prob_join_xy=input;
    quan_size=quan;
    max_run=Max_run;
}

void IB_kernel::smIB()
{

}

 std::vector<std::vector<double>> quantize_to_xt(std::vector<std::vector<double>> & input, std::vector<unsigned> & cluster)
 {
     std::vector<std::vector<double>> join_prob_xt(2,std::vector<double>(cluster.size(),0));
     std::vector<double>::iterator iter1=input[0].begin();
     std::vector<double>::iterator iter2=input[1].begin();
     int i=0;
     for (const unsigned & cluval: cluster)
     {
         join_prob_xt[0][i]=std::accumulate(iter1,iter1+cluval,0);
         join_prob_xt[1][i]=std::accumulate(iter2,iter2+cluval,0);
         iter1=iter1+cluval;
         iter2=iter2+cluval;
     }
     return join_prob_xt;
 }

std::vector<unsigned> IB_kernel::random_cluster()
{
     //initialize random seed
    std::random_device rd;
    std::mt19937 generator(rd());
    //srand(time(0));
    unsigned total_num=prob_join_xy[0].size() ;
    //Para Initial
    std::vector<unsigned> remaining(prob_join_xy[0].size() - 1);
    for (unsigned ii = 0; ii < total_num - 1; ii++)
    {
        remaining[ii] = ii + 1;
    }

    for (const auto &term : remaining)
        std::cout << term;
    std::cout << std::endl;
    std::vector<unsigned> chosen(quan_size - 1);
    int selected;

    //start
    for (unsigned ii = 0; ii < quan_size - 1; ii++)
    {
        std::uniform_int_distribution<> distribution(0,total_num - 2 - ii);
        selected = distribution(generator);
        std::cout << "round: " << ii << "select:" << selected;
        chosen[ii] = remaining[selected];
        std::cout << "number" << chosen[ii] << std::endl;
        remaining.erase(remaining.begin() + selected);
    }

    chosen.push_back(0);
    chosen.push_back(total_num);
    std::sort(chosen.begin(), chosen.end());

    std::vector<unsigned> cluster(quan_size);
    for (unsigned ii = 0; ii < quan_size; ii++)
    {
        cluster[ii] = chosen[ii + 1] - chosen[ii];
    }

    return cluster;
}

