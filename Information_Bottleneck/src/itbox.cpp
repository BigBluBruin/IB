/* information theory box */
#include "Information_Bottleneck/itbox.h"

void flow_stop(double & in)
{
    if(in<pow(10,-15))
        in=pow(10,-15);
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

void ave_joinprob(std::vector<std::vector<double>> & joinprob)
{
    double sum=0;
    std::vector<std::vector<double>>::iterator ounner;
    std::vector<double>::iterator inner;
    for(ounner=joinprob.begin();ounner!=joinprob.end();ounner++)
    {
        for(inner=(*ounner).begin();inner!=(*ounner).begin();inner++)
        {
            flow_stop(*inner);
        }
    }
    for(const auto& out: joinprob)
    {
        for(const auto& term: out)
        {
            sum+=term;
        }
    }
        for(ounner=joinprob.begin();ounner!=joinprob.end();ounner++)
    {
        for(inner=(*ounner).begin();inner!=(*ounner).begin();inner++)
        {
            (*inner)=(*inner)/sum;
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
    double result = 0;
    for (const auto &inner : joindist)
    {
        for (const auto &term : inner)
        {
            result += log2(term);
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
        result = (*iter1) * log2((*iter1) / (*iter2));
    }
    return -result;
}

double it_js(double p1, std::vector<double> dist1, double p2, std::vector<double> dist2)
{
    std::vector<double> ave=p1*dist1+p2*dist2; //realize vector average
    ave_prob(ave);
    return p1*it_kl(dist1,ave)+p2*it_kl(dist2,ave);

}