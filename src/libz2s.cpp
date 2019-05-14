/*
 *  Copyright Peter G. Jensen, all rights reserved.
 */

#include "libz2s.h"
#include "ZonotopStrategy.h"

#include <fstream>
#include <string>
#include <iostream>

void* parse_strategy(const char* file)
{
    if(file == nullptr)
    {
        std::cerr << "Input file must be specified" << std::endl;
        exit(-1);
        return nullptr;
    }
    std::string fn(file);
    std::ifstream instream(fn);
    if(instream.fail())
    {
        std::cerr << "Could not open infile for reading : " << fn << std::endl;
        exit(-1);
        return nullptr;
    }
    std::cerr << "PARSING " << fn << std::endl;
    auto zs = new ZonotopStrategy(ZonotopStrategy::parse(instream));
    std::cerr << "DONE PARSING " << fn << std::endl;
    instream.close();
    return zs;
}
void destroy_strategy(void* ptr)
{
    std::cerr << "destroy " << ptr << std::endl;
    delete (ZonotopStrategy*)ptr;
}

int num_patterns(void* ptr)
{
    std::cerr << "num_patterns " << ptr << std::endl;    
    auto pl = ((ZonotopStrategy*)ptr)->num_patterns();
    std::cerr << pl << std::endl;
    return pl;
}

int max_pattern_length(void* ptr)
{
    std::cerr << "max_pattern_length " << ptr << std::endl;
    auto pl = ((ZonotopStrategy*)ptr)->max_pattern_length();
    std::cerr << pl << std::endl;
    return pl;
}

void active(void* ptr, const double* sample, bool* write)
{
//    std::cerr << "active " << ptr << std::endl;
    ((ZonotopStrategy*)ptr)->active(sample, write);
}

int get_pattern(void* ptr, int el, int* write)
{
//    std::cerr << "get_pattern " << ptr << std::endl;
    return ((ZonotopStrategy*)ptr)->get_pattern(el, write);
}

double get_min(void* ptr, int dimen)
{
    return ((ZonotopStrategy*)ptr)->get_min(dimen);
}

double get_max(void* ptr, int dimen)
{
    return ((ZonotopStrategy*)ptr)->get_max(dimen);
}