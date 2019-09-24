
#include "Information_Bottleneck/itbox.h"
#include "Information_Bottleneck/overloadvec.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <random>

/*std::vector<unsigned> random_cluster(const unsigned total_num, const unsigned quan_size)
{
    //initialize random seed
    std::random_device rd;
    std::mt19937 generator(rd());
    //srand(time(0));

    //Para Initial
    std::vector<unsigned> remaining(total_num - 1);
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
}*/

 /*std::vector<std::vector<double>> quantize_to_xt2(std::vector<std::vector<double>> & input, std::vector<unsigned> & cluster)
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
         i+=1;
     }
     return join_prob_xt;
 }*/

int main()
{

    //function testing
    /*std::vector<double> prob(4,-1);
    std::vector<double>::iterator iter;
    for (iter=prob.begin();iter!=prob.end();iter++)
        *iter=1;
    ave_prob(prob);
    for(const double &term: prob)
        std::cout<< term <<"  ";
    std::cout<<std::endl;
    std::cout<<"Entropy: "<<it_entropy(prob);
    std::cout<<std::endl;
    std::vector<std::vector<double>> join(2,std::vector<double>(4,1));
    std::cout<<"before average";
    for (const auto & out: join)
        for (const auto& term: out)
            std::cout<<term<<"  ";
    std::cout<<std::endl;*/

    //operator + testx
    /*std::vector<std::vector<double>> a1{{1,2,3},{4,5,6}};
    std::vector<std::vector<double>> a2{{7,8,9},{2,3,4}};
    std::vector<std::vector<double>> sum;
    std::ofstream myfile("testing_log_jp_sum.txt");
    myfile<<"a1--------------"<<std::endl;
    for(const auto & loo1: a1)
    {
        for (const auto & loo2: loo1)
            myfile<<loo2<<"  ";
        myfile<<std::endl;
    }
    myfile<<"a2--------------"<<std::endl;
    for(const auto & loo1: a2)
    {
        for (const auto & loo2: loo1)
            myfile<<loo2<<"  ";
        myfile<<std::endl;
    }
    myfile<<"exp. result-----"<<std::endl;
    sum=a1+a2;
    for(const auto & loo1: sum)
    {
        for (const auto & loo2: loo1)
            myfile<<loo2<<"  ";
        myfile<<std::endl;
    }
    myfile.close();*/

    //test for randome cluster//
    /*std::ofstream myfile("random_cluster_test.txt");
    myfile << "--------repeat test-------" << std::endl;
    std::vector<unsigned> out = random_cluster(4, 4);
    for (const auto &term : out)
        myfile << term << "  ";
    myfile << std::endl;
    int max_run = 100;
    for (int ii = 0; ii < max_run; ii++)
    {
        out = random_cluster(20, 4);
        myfile << "--------test: " << ii << "-----" << std::endl;
        for (const auto &term : out)
            myfile << term << "  ";
        myfile << std::endl;
    }
    myfile.close();*/
    /*std::vector<std::vector<double>> join{{1,2,3,4,5},{6,7,8,9,10}};
    std::ofstream myfile("quantization_cluster_test.txt");
    myfile<<"original prob-------"<<std::endl;
    for(const auto & firlay: join)
    {
        for(const auto& seclay:firlay)
            myfile<<seclay<<"  ";
        myfile<<std::endl;
    }
    std::vector<unsigned> cluster{1,2,2};
    std::vector<std::vector<double>> quaned=quantize_to_xt2(join,cluster);
    myfile<<"after cluster--------"<<std::endl;
    for(const auto & firlay: quaned)
    {
        for(const auto& seclay:firlay)
            myfile<<seclay<<"  ";
        myfile<<std::endl;
    }*/
    return 0;
}