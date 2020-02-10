/*
 * Copyright (C) 2020 Peter G. Jensen <root@petergjoel.dk>
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
	auto val = ((ZonotopStrategy*)ptr)->get_min(dimen);
std::cerr << "min for " << dimen << " is " << val << std::endl;
    return val;
}

double get_max(void* ptr, int dimen)
{
	auto val = ((ZonotopStrategy*)ptr)->get_max(dimen);
	std::cerr << "max for " << dimen << " is " << val << std::endl;
	return val;
}


void* parse_learned(const char* file)
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
    std::vector<double> exactness;
    auto res = new SimpleTree(SimpleTree::parse(instream, true, true, 0, exactness));
    return res;
}

void destroy_learned(void* ptr)
{
    delete (SimpleTree*)ptr;
}

double weight(void* ptr, const double* disc, const double* cont, int a)
{
    return ((SimpleTree*)ptr)->value(disc, cont, a);
}