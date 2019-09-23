
/*#include "Information_Bottleneck/overloadvec.h"

template <typename T>
std::vector<T> std::vector<T>::operator + (const std::vector<T>& a, const std::vector<T>& b)
{


    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <class T, class Q>
std::vector <T> std::vector<T>::operator * (const Q c, std::vector <T> A)
{
    std::transform (A.begin (), A.end (), A.begin (),
                 std::bind(std::multiplies <T> () , c)) ;
    return A ;
}*/