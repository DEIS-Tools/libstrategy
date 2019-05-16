/* 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *  Copyright Peter G. Jensen, all rights reserved.
 */

/* 
 * File:   SimpleTree.cpp
 * Author: Peter G. Jensen <root@petergjoel.dk>
 * 
 * Created on May 9, 2019, 10:21 PM
 */

#include "SimpleTree.h"
#include "errors.h"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <unordered_set>

using json = nlohmann::json;

SimpleTree SimpleTree::parse(std::istream& input) 
{
    auto raw = json::parse(input);
    if(!raw.is_object())
        throw base_error("Input JSON not well formatted");

    if(!raw["version"].is_number_float() || raw["version"].get<double>() != 1.0)
        throw base_error("Input version not supported");

    if(!raw["type"].is_string() || raw["type"].get<std::string>() != "state->regressor")
        throw base_error("Input type not supported");

    if(!raw["representation"].is_string() || raw["representation"].get<std::string>() != "map")
        throw base_error("Input representation not supported");
    
    if(!raw["actions"].is_object())
        throw base_error("Expected actions object");
    // store actions
    SimpleTree tree;
    tree._actions.resize(raw["actions"].size());
    for(auto it = raw["actions"].begin(); it != raw["actions"].end(); ++it)
    {
        auto id = atoi(it.key().c_str());
        tree._actions[id] = it.value().get<std::string>();
    }
    
    if(!raw["statevars"].is_array())
        throw base_error("Expected statevars as array");
    tree._statevars = raw["statevars"].get<std::vector<std::string>>();
    
    if(!raw["pointvars"].is_array())
        throw base_error("Expected pointvars as array");
    tree._pointvars = raw["pointvars"].get<std::vector<std::string>>();

    
    if(!raw["regressors"].is_object())
        throw base_error("Expected regressors as object");
    auto& regressors = raw["regressors"];
    int minim = -1;
    for(auto it = regressors.begin(); it != regressors.end(); ++it)
    {
        auto key = parse_key(it.key());
        if(key.size() != tree._statevars.size())
        {
            std::cerr << it.key() << std::endl;
            std::cerr << key.size() << std::endl;
            std::cerr << tree._statevars.size() << std::endl;
            throw base_error("Expected cardinality of key does not match cardinality of statevars");
        }
        

        if(!it.value().is_object())
            throw base_error("Expected regressor-element as object");
        auto& reg = it.value();
        if(reg["type"].get<std::string>() != "act->point->val")
            throw base_error("Expected regressor-element of type act->point->val");
        if(reg["representation"].get<std::string>() != "simpletree")
            throw base_error("Expected regressor-element represented as simpletree");
        bool is_minimize = reg["minimize"].is_boolean() ? reg["minimize"].get<bool>() : reg["minimize"].get<int32_t>();
        if(!reg["regressor"].is_object())
            throw base_error("Expected regressor-sub-element as an object");
        if(minim == -1)
            minim = is_minimize;
        if(minim != is_minimize)
            throw base_error("Expected all sub-regressors to have same minimization flag");
        for(auto iit = reg["regressor"].begin(); iit != reg["regressor"].end(); ++iit)
        {
            size_t action = atoi(iit.key().c_str());
            if(action >= tree._actions.size())
                throw base_error("Action-index is out of bounds");
            
            if(tree._root == nullptr)
            {
                tree._root = std::make_shared<node_t>();
                tree._root->_limit = -std::numeric_limits<double>::infinity();
                tree._root->_cost = std::numeric_limits<double>::infinity();
                if(!is_minimize) tree._root->_cost *= -1;
                tree._root->_var = std::numeric_limits<uint32_t>::max();
            }
            auto obj = iit.value();
            tree._root->insert(key, obj, action, tree, 0, is_minimize);
            //tree._root = tree._root->simplify(false, tree._nodemap);
        }
    }
    tree._root->subsumption_reduction(minim, tree);
    tree._root = tree._root->simplify(true, tree._nodemap);
    return tree;
}

std::vector<double> SimpleTree::parse_key(const std::string& key) {
    size_t i = 0;
    std::vector<double> res;
    if(key[i] != '(')
    {
        throw base_error("Key is incorrectly formatted");
    }
    ++i;
    while(i < key.size() && key[i] != ')')
    {
        std::string ok(key.begin() + i, key.end());
        std::istringstream ss(ok);
        {
            std::string num;
            std::getline(ss, num, ',');
            std::istringstream ns(num);
            res.emplace_back();
            ns >> res.back();
            i += num.length() + 1;
        }
    }
    return res;
}

std::pair<double, double> SimpleTree::node_t::compute_min_max() {
    if(_low)
    {
        auto r = _low->compute_min_max();
        _mincost = std::min(r.first, _mincost);
        _maxcost = std::max(r.second, _maxcost);
    }
    if(_high)
    {
        auto r = _high->compute_min_max();
        _mincost = std::min(r.first, _mincost);
        _maxcost = std::max(r.second, _maxcost);        
    }
    if(!std::isinf(_cost) && !std::isnan(_cost))
    {
        _mincost = std::min(_cost, _mincost);
        _maxcost = std::max(_cost, _maxcost);                
    }
    return std::make_pair(_mincost, _maxcost);
}

void SimpleTree::node_t::action_nodes(std::vector<std::shared_ptr<node_t>>& nodes, uint32_t low, uint32_t high, uint32_t varid) {
    if(_var != varid)
    {
        if(!is_leaf() || !(std::isnan(_cost) || std::isinf(_cost)))
            nodes[high] = shared_from_this(); // just to make sure we dont end up working in garbage memory
    }
    else
    {
        if(_low)
            _low->action_nodes(nodes, low, _limit, varid);
        if(_high)
            _high->action_nodes(nodes, _limit+1, high, varid);
    }
}

void SimpleTree::node_t::subsumption_reduction(bool minimization, SimpleTree& parent) 
{
    if(_var < parent._statevars.size())
    {
        if(_low)
            _low->subsumption_reduction(minimization, parent);
        if(_high)
            _high->subsumption_reduction(minimization, parent);
    }
    else 
    {
        if(_var != parent._statevars.size())
        {
            print(std::cerr, 0);
            return;
        }
        std::vector<std::shared_ptr<node_t>> nodes(parent._actions.size());        
        action_nodes(nodes, 0, parent._actions.size()-1, parent._statevars.size());
        auto val = std::numeric_limits<double>::infinity();
        if(!minimization) val *= -1;
        for(auto& n : nodes)
        {
            if(n == nullptr) continue;
            auto mm = n->compute_min_max();
            if(minimization) val = std::min(val, mm.second);
            else             val = std::max(val, mm.first);
        }
        std::vector<std::pair<double,double>> bounds(parent._pointvars.size(), std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()));
        for(auto& n : nodes)
        {
            bool skip = false;
            if(n == nullptr) continue;
            if(minimization && val < n->_mincost)
                skip = true;
            else if(!minimization && val > n->_maxcost)
                skip = true;
            if(skip)
            {
                n->_low = nullptr;
                n->_high = nullptr;
                n->_var = std::numeric_limits<uint32_t>::max();
                n->_limit = std::numeric_limits<double>::infinity();
                n->_cost = std::numeric_limits<double>::infinity();
                if(minimization)
                    n->_cost *= -1; 
                n->_mincost = std::numeric_limits<double>::infinity();
                n->_maxcost = -std::numeric_limits<double>::infinity();
            }
            else
            {
                n->check_tiles(n.get(), nodes, bounds, val, minimization, parent._statevars.size() + 1);
            }
        }
        
        std::set<std::pair<double,node_t*>> values;
        for(auto& n : nodes)
        {
            if(n == nullptr) continue;
            get_ranks(values, n.get());
        }
        std::unordered_map<double,double> replace;
        for(auto& e : values)
        {
            if(replace.count(e.first) > 0) continue;
            double id = replace.size();
            replace[e.first] = id;
//            std::cerr << e.first << " : " << e.second << std::endl; 
        }
        set_ranks(replace);
    }
}

void SimpleTree::node_t::set_ranks(std::unordered_map<double, double>& values) {
    if(is_leaf())
    {
        if(std::isinf(_cost) || std::isnan(_cost)) return;
        _cost = values[_cost];
    }
    else
    {
        if(_low) _low->set_ranks(values);
        if(_high) _high->set_ranks(values);
    }
}


void SimpleTree::node_t::get_ranks(std::set<std::pair<double, node_t*> >& values, node_t* start) 
{
    if(is_leaf())
    {
        if(std::isinf(_cost) || std::isnan(_cost)) return;
        values.emplace(_cost, start);
    }
    else
    {
        if(_low) _low->get_ranks(values, start);
        if(_high) _high->get_ranks(values, start);
    }
}


bool SimpleTree::node_t::subsumes(std::vector<std::pair<double, double> >& bounds, double val, bool minimization, size_t offset) {
    if(is_leaf())
    {
        if(val == _cost) return true;
        if(minimization && val < _cost) return true;
        else if(!minimization && val > _cost) return true;
        return false;
    }
    else
    {
        if(minimization && _maxcost <= val) return true;
        else if(!minimization && _mincost >= val) return true;
        if(bounds[_var-offset].first <= _limit && _low)
            if(!_low->subsumes(bounds, val, minimization, offset))
                return false;
        if(bounds[_var-offset].second > _limit && _high)
            if(!_high->subsumes(bounds, val, minimization, offset))
                return false;
        return true;
    }
}

void SimpleTree::node_t::check_tiles(node_t* start, std::vector<std::shared_ptr<node_t> >& nodes, std::vector<std::pair<double, double> >& bounds, double val, bool minimization, size_t offset) 
{
    bool remove = false;
    if(minimization && val < _mincost)
        remove = true;
    else if(!minimization && val > _maxcost)
        remove = true;

    if(!remove && is_leaf())
    {
        if(!std::isinf(_cost))
        {
            if(minimization && val < _cost)
                remove = true;
            else if(!minimization && val > _cost)
                remove = true;
            else
            {
                for(auto& n : nodes)
                {
                    if(n.get() == start || n == nullptr) continue;
                    if(n->subsumes(bounds, _cost, minimization, offset))
                    {
                        remove = true;
                        break;
                    }
                }
            }
        }
        else remove = true;
    }
    
    if(remove)
    {
        _low = nullptr;
        _high = nullptr;
        _limit = std::numeric_limits<double>::infinity();
        _cost = std::numeric_limits<double>::infinity();
        if(!minimization) _cost *= -1;
        _var = std::numeric_limits<uint32_t>::max();
        compute_min_max();
        return;
    }
    
    if(!is_leaf()) {
        double org = _limit;
        auto& bnd = bounds[_var- offset];
        std::swap(org, bnd.second);
        if(_low)
            _low->check_tiles(start, nodes, bounds, val, minimization, offset);
        std::swap(org, bnd.second);
        std::swap(org, bnd.first);
        if(_high)
            _high->check_tiles(start, nodes, bounds, val, minimization, offset);
        std::swap(org, bnd.first);
        
        if(_low && _low->is_leaf() && _high && _high->is_leaf() && std::isinf(_low->_cost) && std::isinf(_high->_cost))
        {
            _low = nullptr;
            _high = nullptr;
            _cost = std::numeric_limits<double>::infinity();
            if(!minimization) _cost *= -1;
        }
    }
}





void SimpleTree::node_t::insert(std::vector<double>& key, json& tree, size_t action, SimpleTree& parent, size_t prefix, bool minimize)
{    
    if(tree.is_number())
    {
        if(is_leaf())
        {
            if(minimize)
                _cost = std::isinf(_cost) ? tree.get<double>() : std::min(_cost, tree.get<double>());
            else
                _cost = std::isinf(_cost) ? tree.get<double>() : std::max(_cost, tree.get<double>());
        }
        else
        {
            assert(false);
        }
    }
    else if(is_leaf())
    {
        if(std::isinf(_limit))
        {
            assert(std::isinf(_cost));
            assert(_var == std::numeric_limits<uint32_t>::max());
            _low = std::make_shared<node_t>();
            _high = std::make_shared<node_t>();
            if(!minimize)
            {
                _low->_cost *= -1;
                _high->_cost *= -1;
            }
            _low->_parent = _high->_parent = this;
            _var = prefix <= key.size() ? prefix : key.size() + 1 + tree["var"].get<uint32_t>();
            if(prefix == key.size())
                _limit = action;
            else if(prefix < key.size())
                _limit = key[prefix];
            else
                _limit = tree["bound"].get<double>();
            if(prefix <= key.size())
                _low->insert(key, tree, action, parent, prefix+1, minimize);
            else
            {
                _low->insert(key, tree["low"], action, parent, prefix+1, minimize);
                _high->insert(key, tree["high"], action, parent, prefix+1, minimize);
            }
        }
        else
        {
            assert(!std::isinf(_cost));
            assert(_var != std::numeric_limits<uint32_t>::max());
            assert(false);
        }
    }
    else if(prefix == key.size())
    {
        if(_var != prefix)
        {
            //_parent->replace(this, key, tree, action, parent, prefix);
            auto n1 = std::make_shared<node_t>();
            auto n2 = std::make_shared<node_t>();
            if(!minimize)
            {
                n1->_cost *= -1;
                n2->_cost *= -1;
            }
            auto next = n2;
            n2->_parent = n1.get();
            n1->_limit = action;
            n1->_var = prefix;
            n1->_low = n2;
            n1->_high = shared_from_this();
            if(action != 0)
            {
                n2->_limit = action - 1;
                n2->_var = prefix;
                n2->_low = shared_from_this();
                next = n2->_high = std::make_shared<node_t>();
                if(!minimize) next->_cost *= -1;
            }
            if(this == parent._root.get())
            {
                parent._root = n1;
            }
            else
            {
                if(_parent->_low.get() == this) _parent->_low = n1;
                else                            _parent->_high = n1;
                n1->_parent = _parent;
            }
            _parent = n1.get();
            next->insert(key, tree, action, parent, prefix, minimize);
        }
        else
        {
            if(_limit == action)
            {
                _low->insert(key, tree, action, parent, prefix + 1, minimize);
            }
            else
            {
                if(_limit > action)
                    _low->insert(key, tree, action, parent, prefix, minimize);
                else
                    _high->insert(key, tree, action, parent, prefix, minimize);                
            }
        }
    }
    else if(prefix < key.size())
    {
        if(_var != prefix)
        {
            //_parent->replace(this, key, tree, action, parent, prefix);
            assert(false);
        }
        else
        {
            if(_limit == key[prefix])
            {
                _low->insert(key, tree, action, parent, prefix + 1, minimize);
            }
            else
            {
                if(_limit > key[prefix])
                    _low->insert(key, tree, action, parent, prefix, minimize);
                else
                    _high->insert(key, tree, action, parent, prefix, minimize);
            }
        }
    }
    else
    {
        size_t var = tree["var"].get<uint32_t>() + 1 + key.size();
        if(_var != var)
        {
            //_parent->replace(this, key, tree, action, parent, prefix);
            assert(false);
        }
        else
        {
            double bound = tree["bound"].get<double>();
            if(_limit == bound)
            {
                _low->insert(key, tree["low"], action, parent, prefix+1, minimize);
                _high->insert(key, tree["high"], action, parent, prefix+1, minimize);
            }
            else
            {
                if(_limit > bound)
                    _low->insert(key, tree, action, parent, prefix, minimize);
                else
                    _high->insert(key, tree, action, parent, prefix, minimize);
            }
        }
    }
}

double SimpleTree::node_t::value(double* disc, double* cont, uint32_t action) const {
    return _root->value(disc, cont, action, _statevars.size());
}

double SimpleTree::node_t::value(double* disc, double* cont, uint32_t action, size_t ndisc) const {
    if(is_leaf())
        return _cost;
    if(_var < disc)
    {
        if(disc[_var] <= _limit)
            return _low ? _low->value(disc, cont, action, ndisc) : _cost;
        else 
            return _high ? _high->value(disc, cont, action, ndisc) : _cost;
    }
    else if(_var == disc)
    {
        if(action <= _limit)
            return _low ? _low->value(disc, cont, action, ndisc) : _cost;
        else
            return _high ? _high->value(disc, cont, action, ndisc) : _cost;
    }
    else
    {
        if(cont[_var-(disc+1)] <= _limit)
            return _low ? _low->value(disc, cont, action, ndisc) : _cost;
        else
            return _high ? _high->value(disc, cont, action, ndisc) : _cost;        
    }
}



std::shared_ptr<SimpleTree::node_t> SimpleTree::node_t::simplify(bool make_dd, ptrie::map<std::shared_ptr<node_t>>& nodemap) {
    if(_low) _low = _low->simplify(make_dd, nodemap);
    if(_high) _high = _high->simplify(make_dd, nodemap);
    if(std::isinf(_limit) && (_low != nullptr || _high != nullptr))
    {
        return _low == nullptr ? _high : _low;
    }
    
    assert(_low.get() != this);
    assert(_high.get() != this);
    if( (_low == nullptr || _low->is_leaf()) &&
        (_high == nullptr || _high->is_leaf()))
    {
        auto lc = _low ? _low->_cost : _cost;
        auto hc = _high ? _high->_cost : _cost;
        if(lc == hc)
        {
            _low = nullptr;
            _high = nullptr;
            _cost = lc;
            _limit = std::numeric_limits<double>::infinity();
            _var = std::numeric_limits<typeof(_var)>::max();
            return shared_from_this();
        }
    }
    if(make_dd)
    {
        if(_high == nullptr)
            return _low;
        if(_low == nullptr)
            return _high;
        if(_low && _parent && !_low->is_leaf() && _high->is_leaf())
        {
            auto cost = _high ? _high->_cost : _cost;
            auto lcost = _low->_high ? _low->_high->_cost : _low->_cost;
            if(cost == lcost && _low->_var == _var)
                return _low;
        }
        if(_high && _parent && !_high->is_leaf() && _low->is_leaf())
        {
            auto cost = _low ? _low->_cost : _cost;
            auto hcost = _high->_low ? _high->_low->_cost : _high->_cost;
            if(cost == hcost && _high->_var == _var)
                return _high;
        }
    }
    if(_high != nullptr && _high == _low)
    {
        return _high;
    }
    if(make_dd)
    {
        auto sig = signature();
        auto r = nodemap.insert(sig.first.get(), sig.second);
        auto& ptr = nodemap.get_data(r.second);               
        if(ptr == nullptr)
            ptr = nodemap.get_data(r.second) = shared_from_this();
        return ptr;
    }
    return shared_from_this();
}

std::pair<std::unique_ptr<unsigned char[]>, size_t> SimpleTree::node_t::signature() const {
    std::pair<std::unique_ptr<unsigned char[]>, size_t> res;
    res.second = sizeof(uint32_t)*1+sizeof(double)+sizeof(size_t)*2;
    res.first = std::make_unique<unsigned char[]>(res.second);
    memset(res.first.get(), 0, res.second);
    unsigned char* next = res.first.get();    
    if(is_leaf())
    {
        ((double*)next)[0] = _cost;
        res.second = sizeof(double);
    }
    else
    {
        ((uint32_t*)next)[0] = _var;
        next += sizeof(uint32_t);
        ((double*)next)[0] = _limit;
        next += sizeof(double);
        ((size_t*)next)[0] = (size_t)_low.get();
        ((size_t*)next)[1] = (size_t)_high.get();        
    }
    return res;
}



bool SimpleTree::node_t::is_leaf() const {
    return _low == nullptr && _high == nullptr;
}

/*void SimpleTree::node_t::rec_insert(double value, std::vector<std::pair<bool,bool>>& handled, std::vector<std::pair<double, double> >& bounds, bool minimize) {
    auto nid = _var;
    assert(nid < bounds.size());
    bool reset_low = false;
    bool reset_high = false;
    bool best = false;
    if( std::isnan(_cost)            ||
        (minimize && value < _cost)  ||
        (!minimize && value > _cost))
    {
        _cost = value;
        best = true;
    }
    if(bounds[_var].second == _limit && !handled[_var].second)
        reset_high = handled[_var].second = true;
    if(bounds[_var].first == _limit && !handled[_var].first)
        reset_low = handled[_var].first = true;    
    if(bounds[_var].second <= _limit)
    {
        if(handled[_var].first && handled[_var].second)
        {
            if(_var == handled.size()-1) {
                if(best)
                {
                    _low = nullptr;
                    _high = nullptr;
                }
                if(reset_low) handled[_var].first = false;
                if(reset_high) handled[_var].second = false;
                return;
            }
            ++nid;
        }
        if(is_leaf() && !best)
        {
            if(reset_low) handled[_var].first = false;
            if(reset_high) handled[_var].second = false;
            return;
        }
        if(_low == nullptr)
        {
            _low = std::make_shared<node_t>();
            _low->_var = nid;
            _low->_parent = this;
            _low->_limit = (!handled[nid].first ? bounds[nid].first : bounds[nid].second);            
        }
        _low->rec_insert(value, handled, bounds, minimize);
        if(reset_low) handled[_var].first = false;
        if(reset_high) handled[_var].second = false;
    }
    else if(bounds[_var].first >= _limit)
    {
        if(handled[_var].first && handled[_var].second)
        {
            if(_var == handled.size()-1) {
                if(best)
                {
                    _low = nullptr;
                    _high = nullptr;
                }
                if(reset_low) handled[_var].first = false;
                if(reset_high) handled[_var].second = false;
                return;
            }
            ++nid;
            assert(nid < bounds.size());
        }
        if(is_leaf() && !best)
        {
            if(reset_low) handled[_var].first = false;
            if(reset_high) handled[_var].second = false;
            return;
        }
        if(_high == nullptr)
        {
            _high = std::make_shared<node_t>();
            _high->_parent = this;
            _high->_var = nid;
            _high->_limit = (!handled[nid].first ? bounds[nid].first : bounds[nid].second);
        }
        _high->rec_insert(value, handled, bounds, minimize);
        if(reset_low) handled[_var].first = false;
        if(reset_high) handled[_var].second = false;
    }
    else
    {
        if(best && handled[_var].first && handled[_var].second && _var == handled.size()-1) {
            _low = nullptr;
            _high = nullptr;
            if(reset_low) handled[_var].first = false;
            if(reset_high) handled[_var].second = false;
        }
        if(is_leaf() && !best)
        {
            if(reset_low) handled[_var].first = false;
            if(reset_high) handled[_var].second = false;
        }
        if(_low == nullptr)
        {
            _low = std::make_shared<node_t>();
            _low->_var = nid;
            _low->_parent = this;
            _low->_limit = (!handled[nid].first ? bounds[nid].first : bounds[nid].second);            
        }
        if(_high == nullptr)
        {
            _high = std::make_shared<node_t>();
            _high->_parent = this;
            _high->_var = nid;
            _high->_limit = (!handled[nid].first ? bounds[nid].first : bounds[nid].second);            
        }

        _low->rec_insert(value, handled, bounds, minimize);
        _high->rec_insert(value, handled, bounds, minimize);
    }
    if(reset_low) handled[_var].first = false;
    if(reset_high) handled[_var].second = false;
}*/

std::ostream& SimpleTree::print(std::ostream& stream)
{
    _root->print(stream, 0);
    return stream;
}

std::ostream& SimpleTree::node_t::print(std::ostream& out, size_t tabs) const {
    for(size_t i = 0; i < tabs; ++i) out << "\t";
    out << "{\"var\":" << _var << ",\"bound\":" << _limit;
    if(!std::isnan(_cost))
    {
        out << ",\"value\":" << _cost;
    }
    if(_low)
    {
        out << ",\n";
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "\"low\":\n";
        _low->print(out, tabs + 1);
    }
    if(_high)
    {
        out << ",\n";
        for(size_t i = 0; i < tabs; ++i) out << "\t";
        out << "\"high\":\n";
        _high->print(out, tabs + 1);
    }
    out << "\n";
    for(size_t i = 0; i < tabs; ++i) out << "\t";
    out << "}";
    return out;
}

std::ostream& SimpleTree::print_c(std::ostream& stream, std::string name, size_t split) 
{
    
    stream << "double " << name << "(unsigned int action, const bool* active, const double* disc, const double* vars)\n{\n";
    stream << "\t// Depth = " << _root->depth() << std::endl;
    stream << "\t// Size = " << _nodemap.size() << std::endl;
    std::unordered_set<node_t*> printed;
    _root->print_c(stream, _statevars.size(), printed, 1);
/*    for(size_t i = 0; i < _nodemap.size(); ++i)
    {
        auto& node = _nodemap.get_data(i);
        stream << "\tl" << node.get() << ":\n";
        node->print_c(stream, _statevars.size(), 2);
    }*/
    stream << "}\n";
    return stream;
}

size_t SimpleTree::node_t::depth() const {
    if(_low && _high == nullptr)
        return 1 + _low->depth();
    if(_high && _low == nullptr)
        return 1 + _high->depth();
    if(_low == nullptr && _high == nullptr)
        return 0;
    return 1 + std::max(_low->depth(), _high->depth());
}


std::ostream& SimpleTree::node_t::print_c(std::ostream& stream, size_t disc, std::unordered_set<node_t*>& printed, size_t tabs) {

    if(printed.count(this) > 0) return stream;
    printed.insert(this);
    for(size_t i = 0; i < tabs; ++i) stream << "\t";
    stream << "l" << this << ":\n";

    if(is_leaf())
    {
        for(size_t i = 0; i < tabs+1; ++i) stream << "\t";
        stream << "return " << _cost << ";" << std::endl;
        return stream;
    }
    bool print_low = false;
    bool print_high = false;
    if(std::isinf(_limit))
    {
        if(_low)
            return _low->print_c(stream, disc, printed, tabs);
        else if(_high)
            return _high->print_c(stream, disc, printed, tabs);
        return stream;
    }
    for(size_t i = 0; i < tabs+1; ++i) stream << "\t";
    if(_var == disc)
    {
        stream << "if(action <= " << _limit << ") ";
    }
    else
    {
        auto type = disc > _var ? "disc" : "vars";
        auto vid = disc > _var ? _var : _var - (disc+1);
        stream << "if(" << type << "[" << vid << "] <= " << _limit << ") ";
    }
    if(_low)
    {
        if(_low->is_leaf() && std::isinf(_low->_limit))
            stream << "return " << _low->_cost;
        else
        {
            stream << "goto l" << _low.get();
            print_low = true;
        }
    }
    else
        stream << "return " << _cost;
    stream << ";\n";
    for(size_t i = 0; i < tabs+1; ++i) stream << "\t";
    stream << "else ";
    if(_high)
    {
        if(_high->is_leaf() && std::isinf(_high->_limit))
            stream << "return " << _high->_cost;
        else
        {
            stream << "goto l" << _high.get();
            print_high = true;
        }
    }
    else
        stream << "return " << _cost;
    stream << ";\n";
    if(print_low) _low->print_c(stream, disc, printed, tabs);
    if(print_high) _high->print_c(stream, disc, printed, tabs);
    return stream;    
}
