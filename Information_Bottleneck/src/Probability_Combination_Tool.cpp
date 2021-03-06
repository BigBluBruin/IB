#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"

std::vector<std::vector<double>> prob_combination(std::vector<std::vector<double>> &first_input, std::vector<std::vector<double>> &second_input, const char oper_type[])
{
    std::vector<std::vector<double>> aa(2);
    unsigned size_first = first_input[0].size();
    unsigned size_second = second_input[0].size();
    std::vector<std::vector<double>> combined_prob(2, std::vector<double>(size_first * size_second, -1.0));
    if (strcmp(oper_type, "vari") == 0 || strcmp(oper_type, "check") == 0)
    {
        for (unsigned ii1 = 0; ii1 < size_first; ii1++)
        {
            for (unsigned ii2 = 0; ii2 < size_second; ii2++)
            {
                if (strcmp(oper_type, "vari") == 0)
                {
                    combined_prob[0][ii1 * size_second + ii2] = first_input[0][ii1] * second_input[0][ii2] / 0.5;
                    combined_prob[1][ii1 * size_second + ii2] = first_input[1][ii1] * second_input[1][ii2] / 0.5;
                    // if (isnan(combined_prob[0][ii1 * size_second + ii2])||isnan( combined_prob[1][ii1 * size_second + ii2]))
                    // {
                    //     std::cout<<"vari -- NAN"<<first_input[0][ii1]<<" "<<second_input[0][ii2]<<"  "<<first_input[1][ii1]<<"  "<<second_input[1][ii2]<<std::endl;;
                    // }
                }
                else if (strcmp(oper_type, "check") == 0)
                {
                    combined_prob[0][ii1 * size_second + ii2] = first_input[0][ii1] * second_input[0][ii2] + first_input[1][ii1] * second_input[1][ii2];
                    combined_prob[1][ii1 * size_second + ii2] = first_input[0][ii1] * second_input[1][ii2] + first_input[1][ii1] * second_input[0][ii2];
                    //  if (isnan(combined_prob[0][ii1 * size_second + ii2])||isnan( combined_prob[1][ii1 * size_second + ii2]))
                    // {
                    //     std::cout<<"check NAN"<<first_input[0][ii1]<<" "<<second_input[0][ii2]<<"  "<<first_input[1][ii1]<<"  "<<second_input[1][ii2]<<std::endl;;
                    // }
                }
                else
                    std::cout << "Wrong Info: invalide command, plz check " << std::endl;
            }
        }
    }
    else if (strcmp(oper_type, "check_min") == 0)
    {
        combined_prob.assign(2, std::vector<double>(size_first, 0.0));
        if(size_first!=size_second)
        {
            std::cout << "min check wrong... two size not equal cardinality" <<"first size :"<<size_first<<"--second size"<<size_second<< std::endl;
        }
        for (unsigned ii1 = 0; ii1 < size_first; ii1++)
        {
            for (unsigned ii2 = 0; ii2 < size_second; ii2++)
            {
                double prob_0 = first_input[0][ii1] * second_input[0][ii2] + first_input[1][ii1] * second_input[1][ii2];
                double prob_1 = first_input[0][ii1] * second_input[1][ii2] + first_input[1][ii1] * second_input[0][ii2];
                if (isnan(prob_0)||isnan(prob_1))
                {
                    prob_1=0;
                    prob_0=0;
                    std::cout<<"here"<<std::endl;
                }
                if ((ii1) >= (ii2) && (ii2 + 1) > size_first / 2)
                {
                    combined_prob[0][ii2] += prob_0;
                    combined_prob[1][ii2] += prob_1;
                }
                else if ((ii2) >= (ii1) && (ii1 + 1) > size_first / 2)
                {
                    combined_prob[0][ii1] += prob_0;
                    combined_prob[1][ii1] += prob_1;
                }
                else if ((ii1) <= (ii2) && (ii2) <= size_first / 2 - 1)
                {
                    combined_prob[0][size_first - 1 - ii2] += prob_0;
                    combined_prob[1][size_first - 1 - ii2] += prob_1;
                }
                else if ((ii2) <= (ii1) && (ii1) <= size_first / 2 - 1)
                {
                    combined_prob[0][size_first - 1 - ii1] += prob_0;
                    combined_prob[1][size_first - 1 - ii1] += prob_1;
                }
                else if ((ii1 + 1) > size_first / 2 && (ii2 + 1) <= size_first / 2)
                {
                    combined_prob[0][std::max(ii2, size_first - 1 - ii1)] += prob_0;
                    combined_prob[1][std::max(ii2, size_first - 1 - ii1)] += prob_1;
                }
                else if ((ii2 + 1) > size_first / 2 && (ii1 + 1) <= size_first / 2)
                {
                    combined_prob[0][std::max(ii1, size_first - 1 - ii2)] += prob_0;
                    combined_prob[1][std::max(ii1, size_first - 1 - ii2)] += prob_1;
                }
                else
                {
                    std::cout<<"Wrong Info: No such case..." <<std::endl;
                }               
            }
        }
    }

    // if (strcmp(oper_type, "check_min") != 0)
    //{
        for (unsigned ii = 0; ii < combined_prob[0].size(); ii++)
        {
            if (!isnan(combined_prob[0][ii]) && !isnan(combined_prob[1][ii]) && !(combined_prob[0][ii] == 0) && !(combined_prob[1][ii] == 0))
            {
                aa[0].push_back(combined_prob[0][ii]);
                aa[1].push_back(combined_prob[1][ii]);
            }
        }

        combined_prob = aa;
    //}

    return combined_prob;
}

std::vector<std::vector<double>> clip_prob(std::vector<std::vector<double>> &input, double threshold)
{
    std::vector<std::vector<double>> out(2);
    out[0].clear();
    out[1].clear();
    for(unsigned index=0;index<input[0].size();index++)
    {
        if(input[0][index]+input[1][index]>threshold)
        {
            out[0].push_back(input[0][index]);
            out[1].push_back(input[1][index]);
        }
    }
    return out;
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
        iter_start = input[1].begin() + counter;
        iter_end = iter_start + partition[index];
        //std::cout<<partition[index]<<"  "<< aa<<std::endl;
        combined_prob[1].push_back(std::accumulate(iter_start, iter_end, 0.0));

    }
    //std::cout << "finished partition" << std::endl;
    return combined_prob;
}

void prob_offset(std::vector<std::vector<double>> &input, double alpha)
{
    double summ;
    double llr;
    for (unsigned ii = 0 ; ii<input.size(); ii++)
    {
        llr = log(input[ii][0]/input[ii][1]);
        llr = alpha*llr;
        summ = input[ii][0]+input[ii][1];
        input[ii][0] = (summ*exp(llr))/(1+exp(llr));
        input[ii][1] = summ/(1+exp(llr));
    }
}




