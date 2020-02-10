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

/* 
 * File:   ZonotopStrategy.cpp
 * Author: Peter G. Jensen
 * 
 * Created on December 13, 2018, 2:54 PM
 */

#include "ZonotopStrategy.h"

#include <nlohmann/json.hpp>

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>



using json = nlohmann::json;

std::vector<ZonotopStrategy::bound_t> ZonotopStrategy::get_bounds(const std::string& zonotop, size_t& vars)
{
    size_t i = 0;
    size_t vid = 0;
    std::vector<bound_t> bounds;
    while(i < zonotop.size() && zonotop[i] != ')')
    {
        if(zonotop[i] == '[')
        {
            double lower = 0, upper = 0;
            // new element
            ++i;
            std::string ok(zonotop.begin() + i, zonotop.end());
            std::istringstream ss(ok);
            {
                std::string num;
                std::getline(ss, num, ',');
                std::istringstream ns(num);
                ns >> lower;
                i += num.length() + 1;
            }
            while(i < zonotop.size() && zonotop[i] == ' ') ++i;
            {
                std::string num;
                std::getline(ss, num, ']');
                std::istringstream ns(num);
                ns >> upper;        
                i += num.length() + 1;
            }
            bounds.emplace_back(lower, upper);
        }
        ++i;
    }
    if(vars == 0)
        vars = vid;
    if(vars != vid)
    {
        throw base_error("Dimensionality of controller differs");
    }
    return bounds;
}

ZonotopStrategy ZonotopStrategy::parse(std::istream& input) 
{
    auto raw = json::parse(input);
    ZonotopStrategy strategy;
    size_t vars = 0;
    if(raw.is_object())
    {
        for(auto it = raw.begin(); it != raw.end(); ++it)
        {
            auto bounds = get_bounds(it.key(), vars);

            if(it.value().is_array())
            {
                for(auto e : it.value())
                {
                    std::vector<uint16_t> pat = e.get<std::vector<uint16_t>>();
                    if(pat.size() > 0)         
                    {
                        auto res = strategy._states.insert((unsigned char*)pat.data(), pat.size()*sizeof(uint16_t));
                        strategy._max_length = std::max(strategy._max_length, pat.size());
                        if(strategy.add(bounds, res.second))
                        {
                            //std::cerr << "COL " << it.key() << " : " << e << std::endl;
                            //exit(-1);
                        }
                    }
                }
            }
            else
            {
                throw base_error("Input JSON not well formatted");
            }            
        }
    }
    else
    {
        throw base_error("Input JSON not well formatted");
    }
    return strategy;
}

bool ZonotopStrategy::add(std::vector<bound_t>& bounds, size_t state) 
{
    if(_root == nullptr)
    {
        _root = std::make_shared<node_t>();
        _root->_varid = 0;
        _root->_limit = bounds[0]._lower;
    }
    std::vector<std::pair<bool,bool>> handled(bounds.size());
    return rec_insert(_root.get(), bounds, handled, state);
}

bool ZonotopStrategy::rec_insert(node_t* node, std::vector<bound_t>& bounds, std::vector<std::pair<bool,bool>>& handled, size_t state) {
    auto nid = node->_varid;
    assert(nid < bounds.size());
    bool reset = false;
    bool col = false;
    if(bounds[node->_varid]._upper <= node->_limit)
    {
        if(bounds[node->_varid]._upper == node->_limit)
            reset = handled[node->_varid].second = true;
        if(handled[node->_varid].first && handled[node->_varid].second)
        {
            if(node->_varid == bounds.size()-1)
            {
                auto lb = std::lower_bound(node->_low_patterns.begin(), node->_low_patterns.end(), state);
                if(lb != node->_low_patterns.end() && *lb == state)
                {
                    col = true;
                }
                else
                {
                    node->_low_patterns.insert(lb, state);
                    try_merge(node, state);
                }
                return col;
            }
            ++nid;
            assert(nid < bounds.size());
        }

        if(node->_low == nullptr)
        {
            node->_low = std::make_shared<node_t>();
            node->_low->_varid = nid;
            node->_low->_parent = node;
            node->_low->_limit = (!handled[nid].first ? bounds[nid]._lower : bounds[nid]._upper);
        }
        col |= rec_insert(node->_low.get(), bounds, handled, state);
        if(reset) handled[node->_varid].second = false;
    }
    else if(bounds[node->_varid]._lower >= node->_limit)
    {
        if(bounds[node->_varid]._lower == node->_limit)
            reset = handled[node->_varid].first = true;
        if(handled[node->_varid].first && handled[node->_varid].second)
        {
            if(node->_varid == bounds.size()-1)
            {
                auto lb = std::lower_bound(node->_high_patterns.begin(), node->_high_patterns.end(), state);
                if(lb != node->_high_patterns.end() && *lb == state)
                {
                    col = true;
                }
                else
                {
                    node->_high_patterns.insert(lb, state);
                    try_merge(node, state);
                }
                return col;
            }
            ++nid;
            assert(nid < bounds.size());
        }
        
        if(node->_high == nullptr)
        {
            node->_high = std::make_shared<node_t>();
            node->_high->_parent = node;
            node->_high->_varid = nid;
            node->_high->_limit = (!handled[nid].first ? bounds[nid]._lower : bounds[nid]._upper);
        }
        col |= rec_insert(node->_high.get(), bounds, handled, state);
        if(reset) handled[node->_varid].first = false;
    }
    else
    {
        if(node->_low == nullptr)
        {
            node->_low = std::make_shared<node_t>();
            node->_low->_varid = nid;
            node->_low->_parent = node;
            node->_low->_limit = (!handled[nid].first ? bounds[nid]._lower : bounds[nid]._upper);            
        }
        if(node->_high == nullptr)
        {
            node->_high = std::make_shared<node_t>();
            node->_high->_parent = node;
            node->_high->_varid = nid;
            node->_high->_limit = (!handled[nid].first ? bounds[nid]._lower : bounds[nid]._upper);            
        }

        col |= rec_insert(node->_low.get(), bounds, handled, state);
        col |= rec_insert(node->_high.get(), bounds, handled, state);
    }
    return col;
}

void ZonotopStrategy::try_merge(node_t* node, size_t state) {
    if(node == nullptr) return;
    auto hlb = std::lower_bound(node->_high_patterns.begin(), node->_high_patterns.end(), state);
    if(hlb != std::end(node->_high_patterns) && *hlb == state)
    {
        auto llb = std::lower_bound(node->_low_patterns.begin(), node->_low_patterns.end(), state);
        if(llb != std::end(node->_low_patterns) && *llb == state)
        {
            node->_high_patterns.erase(hlb);
            node->_low_patterns.erase(llb);
            bool empty = (node->_high == nullptr && node->_low == nullptr && node->_high_patterns.empty() && node->_low_patterns.empty());
            auto parent = node->_parent;
            if(node->_parent->_low.get() == node)
            {
                parent->_low_patterns.push_back(state);               
                if(empty) parent->_low = nullptr;
            }
            else
            {
                parent->_high_patterns.push_back(state);
                if(empty) parent->_high = nullptr;
            }
            try_merge(parent, state);
        }
    }
}


std::ostream& ZonotopStrategy::node_t::print(std::ostream& out, const ZonotopStrategy* parent, size_t tabs) const {
    for(size_t i = 0; i < tabs; ++i) out << "\t";
    out << "{\"var\":" << _varid << ",\"bound\":" << _limit;
    auto buffer = std::make_unique<uint16_t[]>(parent->_max_length);
    if(!_low_patterns.empty())
    {
        out << ",\n";
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "\"low_patterns\":[";
        bool fp = true;
        for(auto p : _low_patterns)
        {
            if(!fp)
                out << ",";
            out << "[";
            /*size_t length = parent->_states.unpack(p, (unsigned char*)buffer.get());
            bool fe = true;
            for(size_t i = 0; i < length/sizeof(uint16_t); ++i)
            {
                if(!fe)
                    out << ",";
                fe = false;
                out << buffer[i];
            }*/
            out << p;
            out << "]";
        }
        out << "]";
    }
    if(!_high_patterns.empty())
    {
        out << ",\n";
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "\"high_patterns\":[";
        bool fp = true;
        for(auto p : _high_patterns)
        {
            if(!fp)
                out << ",";
            out << "[";
            size_t length = const_cast<ZonotopStrategy*>(parent)->_states.unpack(p, (unsigned char*)buffer.get());
            bool fe = true;
            for(size_t i = 0; i < length/sizeof(uint16_t); ++i)
            {
                if(!fe)
                    out << ",";
                fe = false;
                out << buffer[i];
            }
            out << "]";
        }
        out << "]";
    }
    if(_low)
    {
        out << ",\n";
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "\"low\":\n";
        _low->print(out, parent, tabs + 1);
    }
    if(_high)
    {
        out << ",\n";
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "\"high\":\n";
        _high->print(out, parent, tabs + 1);
    }
    out << "\n";
    for(size_t i = 0; i < tabs; ++i) out << "\t";
    out << "}";
    return out;
}

void ZonotopStrategy::active(const double* sampel, bool* write) const {
    //std::cerr << "LOOKUP [" << sampel[0] << ", " << sampel[1] << "]" << std::endl;
    memset(write, 0, sizeof(bool)*_states.size());
    if(_root == nullptr)
        return;
    auto current = _root.get();
    int active = 0;
//    _root->print(std::cerr, (ZonotopStrategy*)this);
    while(current)
    {
        if(sampel[current->_varid] >= current->_limit)
        {
            for(auto p : current->_high_patterns)
            {
                write[p] = true;
                //std::cerr << "ACTIVE " << p << std::endl;
            }
            active += current->_high_patterns.size();
        }
        if(sampel[current->_varid] <= current->_limit)
        {
            for(auto p : current->_low_patterns)
            {
                write[p] = true;
                //std::cerr << "ACTIVE " << p << std::endl;
            }
            active += current->_low_patterns.size();
        }
        if(sampel[current->_varid] <= current->_limit && current->_low)
        {
            current = current->_low.get();
            continue;
        }
        if(sampel[current->_varid] >= current->_limit && current->_high)
        {
            current = current->_high.get();
            continue;
        }        
        break;
    }
    if(active == 0)
        std::cerr << "No active for [" << sampel[0] << ", " << sampel[1] << "]" << std::endl;
}

int ZonotopStrategy::max_pattern_length() const {
    return _max_length;
}

int ZonotopStrategy::num_patterns() const {
   return _states.size();
}

int ZonotopStrategy::get_pattern(int el, int* write) {
    auto buffer = std::make_unique<uint16_t[]>(_max_length);
    size_t length = _states.unpack(el, (unsigned char*)buffer.get());
    for(size_t i = 0; i < length/sizeof(uint16_t); ++i)
    {
        write[i] = buffer[i];
    }
    return length/sizeof(uint16_t);
}

std::ostream& ZonotopStrategy::print(std::ostream& stream) const {
    if(_root)
        _root->print(stream, this);
    return stream;
}

double ZonotopStrategy::get_max(size_t dimen) const {
    return _root->get_max(dimen);
}

double ZonotopStrategy::get_min(size_t dimen) const {
    return _root->get_min(dimen);
}

double ZonotopStrategy::node_t::get_max(size_t dimen) const {
    double maxval = _varid == dimen ? _limit : -std::numeric_limits<double>::infinity();
    if(_high)
        maxval = std::max(maxval, _high->get_max(dimen));
    if(_low)
        maxval = std::max(maxval, _low->get_max(dimen));
    return maxval;
}

double ZonotopStrategy::node_t::get_min(size_t dimen) const {
    double minval = _varid == dimen ? _limit : std::numeric_limits<double>::infinity();
    if(_high)
        minval = std::min(minval, _high->get_min(dimen));
    if(_low)
        minval = std::min(minval, _low->get_min(dimen));
    return minval;
}



std::ostream& ZonotopStrategy::print_c(std::ostream& stream, std::string function_name) const
{
/*    stream << "const int " << function_name << "_patterns[][" << (_max_length + 1) << "] = {\n";
    bool first = true;
    auto buffer = std::make_unique<uint16_t[]>(_max_length);
    for(size_t i = 0; i < _states.size(); ++i)
    {
        memset(buffer.get(), 0, max_pattern_length()*sizeof(uint16_t));
        auto size = _states.unpack(i, (unsigned char*)buffer.get());
        if(!first)
            stream << ",";
        stream << "{";
        bool bf = true;
        for(int j = 0; j < max_pattern_length()+1; ++j)
        {
            if(!bf)
                stream << ",";
            if((size_t)j < size/sizeof(uint16_t))
                stream << buffer[j];
            else
                stream << "-1";
            bf = false;
        }
        stream << "}\n";
        first = false;
    }    
    stream << "};\n\n";*/
    stream << "bool " << function_name << "(const double* args, bool* patterns)\n";
    stream << "{\n";
    _root->print_c(stream, 1);
    stream << "\treturn false;\n";
    stream << "}\n";
    return stream;
}

std::ostream& ZonotopStrategy::node_t::print_c(std::ostream& out, size_t tabs) const {  
    if(_low != nullptr || !_low_patterns.empty())
    {
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "if(args[" << _varid << "] <= " << _limit << ") {\n";
        if(_low_patterns.size() > 0)
        {
            for(size_t i = 0; i < tabs+1; ++i) out << "\t";
            for(auto p : _low_patterns)
                out << "patterns[" << p << "] = true; ";
            out << "\n";
        }
        if(_low)
            _low->print_c(out, tabs + 1);
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "}\n";
    }

    if(_high != nullptr || !_high_patterns.empty())
    {
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "if(args[" << _varid << "] >= " << _limit << ") {\n";
        if(_high_patterns.size() > 0)
        {
            for(size_t i = 0; i < tabs+1; ++i) out << "\t";
            for(auto p : _high_patterns)
                out << "patterns[" << p << "] = true; ";
            out << "\n";
        }
        if(_high)
            _high->print_c(out, tabs + 1);
        for(size_t i = 0; i < tabs; ++i) out << "\t";    
        out << "}\n";
    }
    return out;
}
