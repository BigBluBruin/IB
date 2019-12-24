
#include "Information_Bottleneck/itbox.h"
#include "Information_Bottleneck/overloadvec.h"
#include "Information_Bottleneck/stats.h"
#include "Information_Bottleneck/IB_kernel.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include "Discrete_Density_Evolution/Regular_DE.h"
#include "Discrete_Density_Evolution/Irregular_DE.h"
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

 /*-std::vector<std::vector<double>> quantize_to_xt2(std::vector<std::vector<double>> & input, std::vector<unsigned> & cluster)
 {
     std::vector<std::vector<double>> join_prob_xt(2,std::vector<double>(cluster.size(),0.0));
     std::vector<double>::iterator iter1=input[0].begin();
     std::vector<double>::iterator iter2=input[1].begin();
     int i=0;
     for (const unsigned & cluval: cluster)
     {
         
         join_prob_xt[0][i]=std::accumulate(iter1,iter1+cluval,0.0);
         join_prob_xt[1][i]=std::accumulate(iter2,iter2+cluval,0.0);
         iter1=iter1+cluval;
         iter2=iter2+cluval;
         i+=1;
     }
     std::cout<<std::endl;
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


    //test: KL & JS divergence 
    /*std::vector<double> p1{0.3,0.3,0.4};
    std::vector<double> p2{0.1,0.5,0.4};
    std::ofstream myfile("testing_KL_JS");
    for (const auto & term: p1)
        myfile<<term<<"  ";
    myfile<<std::endl;
    for (const auto & term: p2)
        myfile<<term<<"  ";
    myfile<<std::endl;
    myfile<<"KL divergence between p1 and p2: "<<it_kl(p1,p2)<<std::endl;
    myfile<<"JS divergence between p1 and p2 with para 0.4 and 0.6: "<<it_js(0.4,p1,0.6,p2);*/

    /*this section test gaussian partition*/
    /*std::vector<std::vector<double>> prob_join_xy=gaussian_disretization(-2,2,6,0.2);
    std::ofstream myfile("testing_channel_dist.txt");
    myfile<<"------min=-2, max=2, sigma^2=0.6, cardi=4------"<<std::endl;
    double sum=0;
    for(const auto & layer1: prob_join_xy)
    {
        for(const auto & term: layer1)
        {
            myfile<<term<<"  ";
            sum+=term;
        }
        myfile<<std::endl;
    }
    myfile<<"probablity sum is: "<<sum<<"."<<std::endl;
    myfile.close();*/
    /*std::vector<double> aa{1,2,3,4,5};
    std::cout<<std::accumulate(aa.begin(),aa.begin()+1,0)<<std::endl;
    std::cout<<std::accumulate(aa.begin(),aa.end(),0)<<std::endl;*/

    /*This part is used to test information bottleneck part*/
    /*double sigma2=0;
    std::vector<double> last_llr;
    for (unsigned index = 0; index < 8; index++)
    {
        sigma2+=0.1;
        unsigned quansize = 16;
        std::vector<std::vector<double>> prob_join_xy = gaussian_disretization(-2,2,2000, sigma2);
        //std::cout<<prob_join_xy[0][0]<<"  "<<prob_join_xy[1][0]<<std::endl;
        IB_kernel kernel_instance(prob_join_xy, quansize, 6000);
        //std::cout<<"------sigma2: "<<sigma2<<"----------"<<std::endl;
        if(index==0)
        {
            kernel_instance.smIB();
            last_llr=llr_cal(kernel_instance.prob_join_xt); 
             for(const auto & bb: last_llr)
                std::cout<<bb<<"  ";
            std::cout<<std::endl;
        }
        else
        {
            kernel_instance.smIB2(last_llr);
            last_llr=llr_cal(kernel_instance.prob_join_xt); 
            for(const auto & bb: last_llr)
                std::cout<<bb<<"  ";
            std::cout<<std::endl;
        }
               
        /*for (const auto &aa : kernel_instance.cluster)
            std::cout << aa << "  ";
        std::cout << std::endl;
        std::cout<<kernel_instance.mi<<std::endl;*/
        /*for(const auto & aa:kernel_instance.prob_join_xt)
        {
            for(const auto & bb: aa)
                std::cout<<bb<<"  ";
            std::cout<<std::endl;
        }
    }*/

    //std::vector<unsigned> partit{170,47,438,57,102,128,21,37,37,21,128,102,57,438,47,170};
    /*std::vector<std::vector<double>> prob_join_xt=quantize_to_xt2(prob_join_xy,partit);
    std::cout<<"------first 170 elements----"<<std::endl;
    for (int ii=0 ;ii<170;ii++)
    {
        std::cout<<prob_join_xy[0][ii]<<"  ";
    }
    std::cout<<std::endl;*/

    /*for(const auto& layer1: prob_join_xt)
    {
        for(const auto&term:layer1)
            std::cout<<term<<"  ";
        std::cout<<std::endl;
    }*/


    /*This part is used to test mutual information*/
    /*std::vector<std::vector<double>> joint{{0.15,0.2,0.15},{0.2,0.1,0.2}};
    std::cout<<it_mi(joint)<<std::endl;*/


    /*This part is used to test index sort*/
    /*std::vector<double> test;
    test.push_back(0.2);
    test.push_back(0.3);
    test.push_back(-0.1);
    test.push_back(0.5);
    std::vector<unsigned> soindex=sorted_pos(test);
    std::cout<<"sorted pos should be: 2-1-0-3-"<<std::endl;
    std::cout<<"result is: ";
    for (const auto &a : soindex)
        std::cout<<a<<"-";
    std::cout<<std::endl;*/

    /*This part is used to test LLR combination*/
    /*std::vector<std::vector<double>> first_input;
    first_input.push_back({1,2,3,4});
    first_input.push_back({5,6,7,8});
    std::vector<std::vector<double>> second_input;
    second_input.push_back({9,10,11,12});
    second_input.push_back({13,14,15,16});
    std::vector<std::vector<double>> check_combin=prob_combination(first_input,second_input,"check");
    std::cout<<"-----check output----"<<std::endl;
    for(const auto &outer: check_combin)
    {
         for(const auto & inner: outer)
            std::cout<<inner<<"  ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;


    std::vector<std::vector<double>> vari_combin=prob_combination(first_input,second_input,"vari");
    std::cout<<"-----vari output----"<<std::endl;
    for(const auto &outer: vari_combin)
    {
         for(const auto & inner: outer)
            std::cout<<inner<<"  ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;*/
    

    /*double dc=6;
    double dv=3;
    std::vector<double> check_edge_dist{0,0,0,0,0,0,0.8140,0.1860};
    std::vector<double> vari_edge_dist{0,0.2558,0.3140,0.0465,0,0,0,0,0,0,0.3837};
    double code_rate=0.5;
    double eb_no_low = 1.0;
    double eb_no_high = 1.5;
    double cur_eb_no;
    unsigned int quansize =16;
    double threshold = pow(10.0, -7.0);
    unsigned max_iter = 50;
    bool flag=true;
    int regular_result;
    double llr_intervel=0.5;
    unsigned int ib_runtime=6000;
    while (flag)
    {
        cur_eb_no =1.25;
        std::cout<<"---------------change ebno as : "<<cur_eb_no<<"-------------------"<<std::endl;
        double sigma2 = pow(10, (-0.1 * cur_eb_no) / (2.0 * code_rate));
        Irregular_DE irregular_ins(check_edge_dist, vari_edge_dist,sigma2,max_iter,quansize,threshold,llr_intervel,ib_runtime);
        regular_result = irregular_ins.Discrete_Density_Evolution();
        /*Regular_DE regularde_ins(dc, dv, sigma2, max_iter, quansize, threshold);
        regular_result = regularde_ins.Discrete_Density_Evolution();*/
       /* switch (regular_result)
        {
        case 1:
            eb_no_high = cur_eb_no;
            break;
        case 2:
            eb_no_low = cur_eb_no;
            break;
        case 3:
            flag=false;
        default:
            break;
        }  
        return 0;     
    }*/

    return 0;
}