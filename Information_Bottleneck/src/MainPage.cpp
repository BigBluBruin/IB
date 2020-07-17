
#include "Information_Bottleneck/itbox.h"
#include "Information_Bottleneck/overloadvec.h"
#include "Information_Bottleneck/stats.h"
#include "Information_Bottleneck/IB_kernel.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include "Discrete_Density_Evolution/Regular_DE.h"
#include "Discrete_Density_Evolution/Irregular_DE.h"
#include "Discrete_Density_Evolution/Quantize_Continuous_DE.h"
#include "Discrete_Density_Evolution/ME_PBRL_DE.h"
#include "Discrete_Density_Evolution/ME2_PBRL_DE.h"
#include <fstream>
#include <time.h>
#include <random>



int main(int argc, char* argv[])
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
        IB_kernel kernel_instance(prob_join_xy, quansize, 6000);
        kernel_instance.Progressive_MMI();         
        for (const auto &aa : kernel_instance.cluster)
            std::cout << aa << "  ";
        std::cout << std::endl;
        //std::cout<<kernel_instance.mi<<std::endl;
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
    


    //double puncture_rate=0;
    //argv[1]: design eb_no
    //argv[2]: design eb_no string   

    //-----80211 LDPC code--------------------------
    // double puncture_rate=0;
    // std::vector<double> check_edge_dist{0,0,0,0,0,0,0.8140,0.1860};
    // std::vector<double> vari_edge_dist{0,0.2558,0.3140,0.0465,0,0,0,0,0,0,0.3837};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000}; //fake
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000}; // fake
    // double code_rate=0.5;

 
    

    //------SA rate 8_9 --------
    // double puncture_rate=0;
    // std::vector<double> check_edge_dist{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    // std::vector<double> vari_edge_dist{0,0,0.375000000000000,0.625000000000000};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000}; //fake
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000}; // fake
    // double code_rate=8.0/9.0;


    //-----kasra pbrl ------------
    // double puncture_rate=0.0588; //0.0588   
    // std::vector<double> vari_edge_dist{0.0833333333333333,0,0,0.142857142857143,0.0595238095238095,0.0714285714285714,0,0,0.107142857142857,0.357142857142857,0,0,0,0,0.178571428571429};  
    // std::vector<double> check_edge_dist{0,0,0.0357142857142857,0,0,0,0.333333333333333,0.190476190476190,0,0,0,0,0,0,0,0,0,0.214285714285714,0.226190476190476};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000};
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000};

    //----rate  1/2 taken from min-LUT paper-----
    // double puncture_rate=0;
    // std::vector<double> vari_edge_dist{0,0.240730062426354,0.210594930679169,0.0300849900970241,0.125354125404267,0,0.0175495775565974,0,0,0,0,0,0,0,0.375686313836588};  
    // std::vector<double> check_edge_dist{0,0,0,0,0,0,0.0216,0.9762,0.0022};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000};
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000};

    //-----rate 1/2 MIM-BP paper -----
    // double puncture_rate=0;
    // std::vector<double> vari_edge_dist{0,0.138119246958694,0.401174892455710,0,0,0,0,0,0.0266432304916971,0,0,0,0,0,0,0,0.434062630093899};  
    // std::vector<double> check_edge_dist{0,0,0,0,0,0,0,0.328969887598871,0.667329663721726,0.00370044867940238};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000};
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000};

    //----(7,4) hamming code----
    // double puncture_rate=0;
    // std::vector<double> vari_edge_dist{0.25,0.5,0.25};  
    // std::vector<double> check_edge_dist{0,0,0,1};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000};
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000};

    //----Regular LDPC----------
    //std::vector<double> vari_edge_dist{0,0,0,1};
    //std::vector<double> check_edge_dist{0,0,0,0,0,0,1};

    //------SA PBRL Code ------------
    // double puncture_rate=0;
    // std::vector<double> check_edge_dist{0,0,0,0,0,0.0681818181818182,0.0795454545454545,0.272727272727273,0.102272727272727,0.113636363636364,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.363636363636364};
    // std::vector<double> vari_edge_dist{0.0795454545454545,0,0.0340909090909091,0,0.0568181818181818,0.0681818181818182,0,0.0909090909090909,0.102272727272727,0.113636363636364,0.250000000000000,0,0,0,0,0,0,0.204545454545455};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000}; //fake
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000}; // fake
    // double code_rate=0.5;

    //-----Nake SA PBRL -------------
    // double puncture_rate=0;
    // std::vector<double> check_edge_dist{0,0,0,0.0617283950617284,0.0740740740740741,0.259259259259259,0.0987654320987654,0.111111111111111,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.395061728395062};
    // std::vector<double> vari_edge_dist{0,0,0.0370370370370370,0,0.0617283950617284,0.0740740740740741,0,0.0987654320987654,0.111111111111111,0.123456790123457,0.271604938271605,0,0,0,0,0,0,0.222222222222222};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000}; //fake
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000}; // fake
    // double code_rate=0.5;


     //------SA rate 8_9 --------
    // double puncture_rate=0;
    // std::vector<double> check_edge_dist{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    // std::vector<double> vari_edge_dist{0,0,0,1};
    // std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000}; //fake
    // std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000}; // fake
    // double code_rate=8.0/9.0;

    //------SS 3K --------
    double puncture_rate=0;
    std::vector<double> check_edge_dist{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6246,0.3754};
    std::vector<double> vari_edge_dist{0, 0.0599, 0 , 0.7634,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1767};
    std::vector<double> eff_check_edge_dist{0,0,0,0,0,0.500000000000000,0.250000000000000,0,0,0,0,0,0,0,0,0,0.250000000000000}; //fake
    std::vector<double> eff_vari_edge_dist{0.821400000000000,0,0,0.0476000000000000,0.131000000000000}; // fake
    double code_rate=7.0/8.0;

    double cur_eb_no;
    unsigned int quansize =32;
    double threshold = pow(10.0, -7.0);
    unsigned max_iter = 50;
    bool flag=true;
    int regular_result;
    double llr_intervel=0.01;
    unsigned int ib_runtime=100;
    std::vector<double> eb_no{std::stod(argv[1])};
    std::vector<std::string> suffix{argv[2]};
    for(unsigned index=0;index<eb_no.size();index++)
    {
        cur_eb_no =eb_no[index];
        double sigma2 = pow((double)10, (-0.1 * cur_eb_no)) / (2.0 * code_rate);
        std::cout<<"---------------change ebno as : "<<cur_eb_no<<"|| corresponding sigma2:"<<sigma2<<"""-------------------"<<std::endl;
        Irregular_DE irregular_ins(check_edge_dist, vari_edge_dist,sigma2,max_iter,quansize,threshold,llr_intervel,ib_runtime,suffix[index],eff_check_edge_dist,eff_vari_edge_dist,puncture_rate);
        if(puncture_rate==0)
        {
            regular_result = irregular_ins.Discrete_Density_Evolution();
        }
        else
        {
            regular_result = irregular_ins.Discrete_Density_Evolution_punc();
        }
    }

    //QCDE("Density_Evolution.txt",16);
    

    //This part tests read description of ME_PBRL_Code
    //argv[1] ->eb_no 
    //argv[2] -> correspoded string
    /*unsigned int quansize = 16;
    double threshold = pow(10.0, -10.0);
    unsigned max_iter = 50;
    bool flag = true;
    int regular_result;
    double llr_intervel = 0.01;
    unsigned int ib_runtime = 50;
    double code_rate = 0.5;
    double cur_eb_no;
    std::vector<double> eb_no{std::stod(argv[1])};

    std::vector<std::string> suffix{argv[2]};
    //std::vector<double> eb_no{0.25,0.27,0.29,0.31,0.33,0.35};
    //std::vector<std::string> suffix{"025", "027", "029", "031", "033"};
    for (unsigned ii = 0; ii < eb_no.size(); ii++)
    {
        cur_eb_no = eb_no[ii];
        std::cout << "---------------change ebno as : " << cur_eb_no << "-------------------" << std::endl;
        double sigma2 = pow(10, (-0.1 * cur_eb_no) / (2.0 * code_rate));
        std::cout << "sigma 2: " << sigma2 << std::endl;
        std::string filename = "PBRL_MET_description.txt";
        ME_PBRL_DE me_pbrl_ins(filename, max_iter, quansize, sigma2, threshold, llr_intervel, ib_runtime, suffix[ii]);
        me_pbrl_ins.read_decription();
        me_pbrl_ins.density_evolution();
    }*/

    //This parts implement ME2_PBRL_DE
    // unsigned int quansize = 32;
    // double threshold = pow(10.0, -10.0);
    // unsigned max_iter = 16;
    // bool flag = true;
    // int regular_result;
    // double llr_intervel = 0.01;
    // unsigned int ib_runtime = 200;
    // double code_rate = 0.5;
    // double cur_eb_no;
    // std::vector<double> eb_no{std::stod(argv[1])};
    // std::vector<std::string> suffix{argv[2]};
    // for (unsigned ii = 0; ii < eb_no.size(); ii++)
    // {
    //     cur_eb_no = eb_no[ii];
    //     std::cout << "---------------change ebno as : " << cur_eb_no << "-------------------" << std::endl;
    //     double sigma2 = pow(10, (-0.1 * cur_eb_no) / (2.0 * code_rate));
    //     std::cout << "sigma 2: " << sigma2 << std::endl;
    //     std::string filename = "PBRL_MET2_description.txt";
    //     ME2_PBRL_DE me2_pbrl_ins(filename, max_iter, quansize, sigma2, threshold, llr_intervel, ib_runtime, suffix[ii]);
    //     me2_pbrl_ins.read_decription();
    //     me2_pbrl_ins.density_evolution();
    // }
    // return 0;
} /*






















//------------------------------recycle bin------------------------------------------------------
 /*switch (regular_result)
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
        return 0;     */
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

/*
   ME_PBRL_DE ins_mepbrl("PBRL_MET_description.txt", 6, 6, 6, 6, 6, "0");
    if (ins_mepbrl.read_decription())
    {
        std::vector<std::vector<double>> cur_output;

        cur_output = ins_mepbrl.variable_node_degree;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.check_node_degree;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.vari_edge_deg_1;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.vari_edge_deg_2;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.vari_edge_deg_3;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.check_edge_deg_1;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.check_edge_deg_2;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
        cur_output = ins_mepbrl.check_edge_deg_3;
        std::cout << "--------------------------" << std::endl;
        for (const auto &outer : cur_output)
        {
            for (const auto &inner : outer)
                std::cout << inner << "  ";
            std::cout << std::endl;
        }
    }
 */

  /*      //-------test can be delete------------------------
        for(const auto aa: me2_pbrl_ins.vari_edge_deg_1)
        {
            for(const auto bb:aa)
                std::cout<<bb<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
        for(const auto aa: me2_pbrl_ins.check_edge_deg_1)
        {
            for(const auto bb:aa)
                std::cout<<bb<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
        //----------------------------------------------------*/