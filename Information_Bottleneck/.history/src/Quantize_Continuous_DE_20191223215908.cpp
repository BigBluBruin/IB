#include "Discrete_Density_Evolution/Quantize_Continuous_DE.h"

bool QCDE(const std::string filename, unsigned quansize)
{
    std::ifstream handle_file(filename);
    if (handle_file.is_open())
    {
        unsigned total_points, max_iter;
        handle_file >> total_points >> max_iter;
        std::vector<double> llr_val(total_points);
        for (unsigned index1 = 0; index1 < total_points; index1++)
        {
            handle_file >> llr_val[index1];
        }
        std::vector<std::vector<double>> check_llr(max_iter);
        for (unsigned index1 = 0; index1 < max_iter; index1++)
        {
            check_llr[index1].assign(total_points, -1);
            for (unsigned index2 = 0; index2 < total_points; index2++)
            {
                handle_file >> check_llr[index1][index2];
            }
        }
        std::vector<std::vector<double>> vari_llr(max_iter);
        for (unsigned index1 = 0; index1 < max_iter; index1++)
        {
            vari_llr[index1].assign(total_points, -1);
            for (unsigned index2 = 0; index2 < total_points; index2++)
            {
                handle_file >> vari_llr[index1][index2];
            }
        }
        std::cout<<"finished written"<<std::endl;
        std::vector<std::vector<double>> joint_prob_xy(2), joint_prob_xt(2);
        std::vector<std::vector<double>> check_threshold,vari_threshold;
        std::vector<std::vector<double>> check_recons, vari_recons;
        std::vector<double> temp;
        double threshold = pow(10, -8);
        double cur_llr;
        for (unsigned cur_iter = 0; cur_iter < max_iter; cur_iter++)
        {
            joint_prob_xy[0].clear();
            joint_prob_xy[1].clear();
            for (unsigned ii = 0; ii < total_points / 2; ii++)
            {
                if (check_llr[cur_iter][ii] > threshold)
                {
                     std::cout<<"been here "<<joint_prob_xy[0].size()<<std::endl;
                    joint_prob_xy[0].push_back((1.0 / (1 + exp(-llr_val[ii]))) * check_llr[cur_iter][ii]);
                    joint_prob_xy[1].push_back(((exp(-llr_val[ii]) / (1 + exp(-llr_val[ii]))) * check_llr[cur_iter][ii]));
                }
                temp = joint_prob_xy[0];
                std::reverse(temp.begin(), temp.end());
                std::copy(temp.begin(), temp.end(), std::back_inserter(joint_prob_xy[0]));
                temp = joint_prob_xy[1];
                std::reverse(temp.begin(), temp.end());
                std::copy(temp.begin(), temp.end(), std::back_inserter(joint_prob_xy[1]));
                //std::cout<<"been here "<<joint_prob_xy[0].size()<<std::endl;
                ave_joinprob_llr(joint_prob_xy, pow(10, -80.0));
                IB_kernel this_check_ib(joint_prob_xy, quansize, 6000);
                this_check_ib.smIB();
                check_threshold.push_back(this_check_ib.threshold);
                check_recons.push_back(llr_cal(this_check_ib.prob_join_xt));
            }
            std::cout<<"finished check iter "<<cur_iter<<std::endl;
        }
        for (unsigned cur_iter = 0; cur_iter < max_iter; cur_iter++)
        {
            joint_prob_xy[0].clear();
            joint_prob_xy[1].clear();
            for (unsigned ii = 0; ii < total_points / 2; ii++)
            {
                if (vari_llr[cur_iter][ii] > threshold)
                {
                    std::cout<<"been here "<<joint_prob_xy[0].size()<<std::endl;
                    joint_prob_xy[0].push_back((1.0 / (1 + exp(-llr_val[ii]))) * vari_llr[cur_iter][ii]);
                    joint_prob_xy[1].push_back(((exp(-llr_val[ii]) / (1 + exp(-llr_val[ii]))) * vari_llr[cur_iter][ii]));
                }
                temp = joint_prob_xy[0];
                std::reverse(temp.begin(), temp.end());
                std::copy(temp.begin(), temp.end(), std::back_inserter(joint_prob_xy[0]));
                temp = joint_prob_xy[1];
                std::reverse(temp.begin(), temp.end());
                std::copy(temp.begin(), temp.end(), std::back_inserter(joint_prob_xy[1]));
                ave_joinprob_llr(joint_prob_xy, pow(10, -80.0));
                IB_kernel this_vari_ib(joint_prob_xy, quansize, 6000);
                this_vari_ib.smIB();
                check_threshold.push_back(this_vari_ib.threshold);
                check_recons.push_back(llr_cal(this_vari_ib.prob_join_xt));
            }
        }
        //Part 1: threholds
        std::ofstream outpufile("threshold.txt");
        if (outpufile.is_open())
        {
            outpufile << quansize << " " << max_iter << "  " << std::endl;
            for (const auto &iter : check_threshold)
            {
                for (const auto &term : iter)
                    outpufile << term << "  ";
                outpufile << std::endl;
            }
            outpufile << std::endl
                      << std::endl;
            for (const auto &iter : vari_threshold)
            {
                for (const auto &term : iter)
                    outpufile << term << "  ";
                outpufile << std::endl;
            }
            outpufile.close();
        }

        //Part 2: Reconstruction function
        std::ofstream handle_reconstrction("reconstruction.txt");
        if (handle_reconstrction.is_open())
        {
            handle_reconstrction << quansize << "  " << max_iter << "  " << std::endl;
            for (const auto &iter : check_recons)
            {
                for (const auto &term : iter)
                {
                    handle_reconstrction << term << "  ";
                }
                handle_reconstrction << std::endl;
            }
            handle_reconstrction << std::endl
                                 << std::endl;
            for (const auto &iter : vari_recons)
            {
                for (const auto &term : iter)
                {
                    handle_reconstrction << term << "  ";
                }
                handle_reconstrction << std::endl;
            }
            handle_reconstrction.close();
        }
        else
        {
            std::cout << "Reconstruction failed to wirte ..." << std::endl;
        }
    }
    else
    {
         std::cout<<"Wrong Info: Didn't see the snow"<<std::endl;
    }
}