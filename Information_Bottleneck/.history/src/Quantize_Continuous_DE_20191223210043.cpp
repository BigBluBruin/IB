#include "Quantize_Continuous_DE.h"

bool QCDE(const std::string filename)
{
    std::ofstream handle_file(filename);
    if(handle_file.is_open())
    {
        unsigned total_points,max_iter;
        handle_file>>total_points>>max_iter;
    }
    else
    {
        
    }
    
}