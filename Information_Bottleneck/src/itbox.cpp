/* information theory box */
#include "Information_Bottleneck/itbox.h"

void flow_stop(double & in)
{
    if(in<pow(10,-100))
        in=pow(10,-100);
}

void ave_prob(std::vector <double> & prob)
{
    std::vector<double>::iterator pr;
    double sum=0;
    for(pr=prob.begin();pr!=prob.end();pr++)
        flow_stop(*pr);
    for(const auto& term:prob)
    {
        sum+=term;
    }
    for(pr=prob.begin();pr!=prob.end();pr++)
    {
        *pr=(*pr)/sum;
    }

}

void ave_joinprob(std::vector<std::vector<double>> &joinprob)
{
    double sum = 0;
    std::vector<std::vector<double>>::iterator ounner;
    std::vector<double>::iterator inner;
    for (ounner = joinprob.begin(); ounner != joinprob.end(); ounner++)
    {
        for (inner = (*ounner).begin(); inner != (*ounner).end(); inner++)
        {
            flow_stop(*inner);
        }
    }

    sum=0;
    for (const auto &out : joinprob)
    {
        for (const auto &term : out)
        {
            sum += term;
        }
    }
    for (ounner = joinprob.begin(); ounner != joinprob.end(); ounner++)
    {
        for (inner = (*ounner).begin(); inner != (*ounner).end(); inner++)
        {
            (*inner) = (*inner) / sum;
        }
    }

}

void ave_joinprob_llr(std::vector<std::vector<double>> &joinprob, double threshold)
{
    double new_1,new_2,sum;
    std::vector<std::vector<double>>::iterator ounner;
    std::vector<double>::iterator inner;
    for (unsigned index = 0; index < joinprob[0].size() / 2; index++)
    {
        if (joinprob[0][index] < threshold && joinprob[1][index] < threshold)
        {
            //Preserves LLR information...
            new_1=threshold;
            new_2=threshold/(joinprob[0][index]/joinprob[1][index]);
            joinprob[0][index]=new_1;
            joinprob[1][index]=new_2;
        }
    }
    for (unsigned index = joinprob[0].size() / 2; index < joinprob[0].size(); index++)
    {
        joinprob[0][index]=joinprob[1][joinprob[0].size()-1-index];
        joinprob[1][index]=joinprob[0][joinprob[0].size()-1-index];
    }
    sum = 0;
    for (const auto &out : joinprob)
    {
        for (const auto &term : out)
        {
            sum += term;
        }
    }
    for (ounner = joinprob.begin(); ounner != joinprob.end(); ounner++)
    {
        for (inner = (*ounner).begin(); inner != (*ounner).end(); inner++)
        {
            (*inner) = (*inner) / sum;
        }
    }
}

double it_entropy(std::vector<double> dist)
{
    double it_ent = 0;
    for (const double & sub_dist : dist)
    {
        it_ent += sub_dist*log2(sub_dist);
    }
    return -it_ent;
}

double it_mi(std::vector<std::vector<double>> joindist)
{
    /*  dimension 1 -> y 
        dimension 2 -> x
        To calculate llr faster */
 
    std::vector<double> probx(2,0);
    for(auto const &term: joindist[0])
        probx[0]+=term;
    for(auto const & term: joindist[1])
        probx[1]+=term;
    
    //std::cout<<"p.m.f. of x:"<<std::endl;
    //std::cout<<probx[0]<<"  "<<probx[1]<<std::endl;
    
    std::vector<double> proby(joindist[0].size(),1);
    for(unsigned ii=0; ii< joindist[0].size(); ii++)
    {
        proby[ii]=joindist[0][ii]+joindist[1][ii];
    }       
    /*std::cout<<"p.m.f. of y:"<<std::endl;
    for(const auto& term: proby)
        std::cout<<term<<"  ";
    std::cout<<std::endl;*/

    double result = 0;

    for(unsigned ii=0; ii<2;ii++)
    {
        for(unsigned jj=0;jj<proby.size();jj++)
        {
            result+=joindist[ii][jj]*log2(joindist[ii][jj]/(probx[ii]*proby[jj]));
        }
    }
    
    return result;
}

double it_kl(std::vector<double> trudist, std::vector<double> appdist)
{
    double result = 0;
    std::vector<double>::iterator iter1, iter2;
    for (iter1 = trudist.begin(), iter2 = appdist.begin(); iter1 < trudist.end() && iter2 < appdist.end(); iter1++, iter2++)
    {
        result += (*iter1) * log2((*iter1) / (*iter2));
    }
    return result;
}

double it_js(double p1, std::vector<double> dist1, double p2, std::vector<double> dist2)
{
    std::vector<double> ave=p1*dist1+p2*dist2; //realize vector average
    ave_prob(ave);
    return p1*it_kl(dist1,ave)+p2*it_kl(dist2,ave);

}

std::vector<double> llr_cal(std::vector<std::vector<double>> & joint_prob)
{
    std::vector <double> llr(joint_prob[0].size(),-1);
    for (unsigned index1=0; index1<=joint_prob[0].size(); index1++)
    {
        llr[index1]=log(joint_prob[0][index1]/joint_prob[1][index1]);
    }
    return llr;
}

void adjust_joint_prob(std::vector<std::vector<double>> & joint_prob, double max_allowed_llr, double allocated_prob)
{
    std::vector<double> llr = llr_cal(joint_prob);
    double new_llr;
    for (unsigned ii =0 ;ii <llr.size();ii++)
    {
        if(abs(llr[ii])>max_allowed_llr)
        {
            if(llr[ii]>0)
            {
                new_llr = max_allowed_llr;
            }
            else
            {
                new_llr = -max_allowed_llr;
            }
            joint_prob[0][ii] = allocated_prob*exp(new_llr)/(1+exp(new_llr));
            joint_prob[1][ii] = allocated_prob/(1+exp(new_llr));
        }
    }
    ave_joinprob(joint_prob);
}


 bool llr_verification(std::vector<double> old_llr, std::vector<double> new_llr)
 {
     for (unsigned ii=0; ii<new_llr.size()/2;ii++)
     {
         if(new_llr[ii]>old_llr[ii])
            return false;
     }
     return true;
 }


 bool cluster_verification(std::vector<unsigned> old_cluster, std::vector<unsigned> new_cluster)
 {
     int sum=0;
     for (unsigned ii=0; ii<new_cluster.size()/2;ii++)
     {
         sum+=abs((int)old_cluster[ii]-(int)new_cluster[ii]);
     }

     if(sum <= 2)
     {
         return false;
     }
     return true;
 }