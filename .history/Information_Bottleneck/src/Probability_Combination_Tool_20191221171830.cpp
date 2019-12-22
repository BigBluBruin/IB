#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"

std::vector<std::vector<double>> prob_combination (std::vector<std::vector<double>> & first_input, std::vector<std::vector<double>> & second_input, const char oper_type[])
{
    unsigned size_first=first_input[0].size();
    unsigned size_second=second_input[0].size();
    std::vector<std::vector<double>> combined_prob(2,std::vector<double>(size_first*size_second,-1.0));
    for (unsigned ii1 = 0; ii1 < size_first; ii1++)
    {
        for (unsigned ii2 = 0; ii2 < size_second; ii2++)
        {
            if(oper_type=="vari")
            {
                combined_prob[0][ii1*size_second+ii2]=first_input[0][ii1]*second_input[0][ii2]/0.5;
                combined_prob[1][ii1*size_second+ii2]=first_input[1][ii1]*second_input[1][ii2]/0.5;

            }
            else if (oper_type=="check")
            {
                combined_prob[0][ii1*size_second+ii2]=first_input[0][ii1]*second_input[0][ii2]+first_input[1][ii1]*second_input[1][ii2];
                combined_prob[1][ii1*size_second+ii2]=first_input[0][ii1]*second_input[1][ii2]+first_input[1][ii1]*second_input[0][ii2];
            }
            else
                std::cout<<"Wrong Info: invalide command, plz check "<<std::endl;
        }
    }
    return combined_prob;
}


std::vector<std::vector<double>> prob_combination_v2 (std::vector<std::vector<double>> & first_input, std::vector<std::vector<double>> & second_input, const char oper_type[], double threshold)
{
    unsigned size_first=first_input[0].size();
    unsigned size_second=second_input[0].size();
    std::vector<std::vector<double>> combined_prob(2);
    double a_sum,upper_value,lower_value;
    for (unsigned ii1 = 0; ii1 < size_first; ii1++)
    {
        for (unsigned ii2 = 0; ii2 < size_second; ii2++)
        {
            if(oper_type=="vari")
            {
                upper_value=first_input[0][ii1]*second_input[0][ii2]/0.5;
                lower_value=first_input[1][ii1]*second_input[1][ii2]/0.5;
                if(upper_value+lower_value>threshold)
                {
                    combined_prob[0].push_back(upper_value);
                    combined_prob[1].push_back(lower_value);
                }

            }
            else if (oper_type=="check")
            {
                combined_prob[0][ii1*size_second+ii2]=first_input[0][ii1]*second_input[0][ii2]+first_input[1][ii1]*second_input[1][ii2];
                combined_prob[1][ii1*size_second+ii2]=first_input[0][ii1]*second_input[1][ii2]+first_input[1][ii1]*second_input[0][ii2];
            }
            else
                std::cout<<"Wrong Info: invalide command, plz check "<<std::endl;
        }
    }
    return combined_prob;
}


std::vector<unsigned> sorted_pos(std::vector<double> &input)
{
    std::vector<unsigned> y(input.size());
    std::size_t n(0);
    std::generate(std::begin(y), std::end(y), [&] { return n++; });

    std::sort(std::begin(y),
              std::end(y),
              [&](unsigned i1, unsigned i2) { return input[i1] < input[i2]; });
    return y;
}

void prob_sort(std::vector<std::vector<double>> & input)
{
    std::vector<double> llr=llr_cal(input);
    std::vector<std::vector<double>> temp=input;
    std::vector<unsigned> pos=sorted_pos(llr);
    for (unsigned ii = 0; ii < llr.size(); ii++)
    {
        input[0][ii]=temp[0][pos[ii]];
        input[1][ii]=temp[1][pos[ii]];
    }
}

std::vector<std::vector<double>> llr_combination(std::vector<std::vector<double>> &input, double threshold)
{
    std::vector<double> llr=llr_cal(input);
    /*for(const auto &aa:llr)
    {
        std::cout<<aa<<" ";
    }*/
    //std::cout<<std::endl;
    std::vector<std::vector<double>> combined_prob(2);
    std::vector<unsigned> partition;
    std::vector<double>::iterator iter_start, iter_end;
    unsigned partition_ind = 0;
    unsigned first_llr_ind = 0;
    unsigned counter = 0;
    partition.push_back(1);
    for (unsigned ii = 1; ii < llr.size() / 2; ii++)
    {
        if ((llr[ii] - llr[first_llr_ind]) < threshold)
        {
            partition[partition_ind]++;
        }
        else
        {
            first_llr_ind = ii;
            partition_ind++;
            partition.push_back(1);
        }
    }

    std::vector<unsigned> partition_reverse = partition;
    std::reverse(partition_reverse.begin(), partition_reverse.end());
    std::copy(partition_reverse.begin(), partition_reverse.end(), std::back_inserter(partition));
    //std::cout<<partition.size()<<std::endl;

    for (unsigned index = 0; index < partition.size(); index++)
    {
        counter=0;
        for (unsigned jj = 0; jj < index; jj++)
            counter += partition[jj];
        iter_start = input[0].begin() + counter;
        iter_end = iter_start + partition[index];
        double aa=std::accumulate(iter_start, iter_end, 0.0);
        //std::cout<<partition[index]<<"  "<< aa<<std::endl;
        combined_prob[0].push_back(aa);
        //std::cout<<"here"<<std::endl;
        iter_start = input[1].begin() + counter;
        iter_end = iter_start + partition[index];
        //std::cout<<partition[index]<<"  "<< aa<<std::endl;
        combined_prob[1].push_back(std::accumulate(iter_start, iter_end, 0.0));

    }
    //std::cout << "finished partition" << std::endl;
    return combined_prob;
}








