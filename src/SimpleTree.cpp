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
        
        auto res = tree._regressors.insert((unsigned char*)key.data(), key.size()*sizeof(double));
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
        bool first = true;
        auto& node = tree._regressors.get_data(res.second);
        for(auto it = reg["regressor"].begin(); it != reg["regressor"].end(); ++it)
        {
            size_t action = atoi(it.key().c_str());
            if(action >= tree._actions.size())
                throw base_error("Action-index is out of bounds");
            
            if(first)
            {
                node._limit = -std::numeric_limits<double>::infinity();
                node._cost = std::numeric_limits<double>::quiet_NaN();
            }
            node.merge(it.value(), action, tree._pointvars.size(), is_minimize);
            first = false;
        
        }
        node.simplify(true);
    }
    for(size_t i = 0; i < tree._regressors.size(); ++i)
    {
        tree._regressors.get_data(i).simplify(false);
    }
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

void SimpleTree::node_t::merge(json tree, size_t action, size_t vars, bool minimize) {
    
    std::vector<std::pair<double,double>> bounds;
    bounds.resize(vars, std::pair<double,double>(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()));
    find_and_insert(tree, action, bounds, minimize);
}

void SimpleTree::node_t::find_and_insert(json tree, size_t action, std::vector<std::pair<double, double> >& bounds, bool minimize) {
    if(tree.is_object())
    {
        size_t var = tree["var"].get<size_t>();
        double ov = tree["bound"].get<double>();
        std::swap(ov, bounds[var].second);
        find_and_insert(tree["low"], action, bounds, minimize);
        std::swap(ov, bounds[var].second);
        std::swap(ov, bounds[var].first);
        find_and_insert(tree["high"], action, bounds, minimize);
        std::swap(ov, bounds[var].first);
    }
    else if(tree.is_number_float())
    {
        // do the merge
        insert(tree.get<double>(), action, bounds, minimize);
    }
    else
    {
        throw base_error("Leafs of trees are expected to be doubles");
    }
}

void SimpleTree::node_t::insert(double value, size_t action, std::vector<std::pair<double, double> >& bounds, bool minimize) 
{
    std::vector<std::pair<bool,bool>> handled(bounds.size());
    rec_insert(value, action, handled, bounds, minimize);
}

void SimpleTree::simplify() {
    for(size_t i = 0; i < _regressors.size(); ++i)
        _regressors.get_data(i).simplify(false);
}


void SimpleTree::node_t::simplify(bool cost_consistent) {
    if(_low) _low->simplify(cost_consistent);
    if(_high) _high->simplify(cost_consistent);
    if(_low && _low->is_leaf() && (_low->_cost == _cost || !cost_consistent) && _low->_action == _action) 
        _low = nullptr;
    if(_high && _high->is_leaf() && (_high->_cost == _cost || !cost_consistent) && _high->_action == _action) 
        _high = nullptr;
    if(_low && _parent && !cost_consistent && (_high == nullptr || (!_low->is_leaf() && _high->is_leaf())))
    {
        auto action = _high ? _high->_action : _action;
        auto laction = _low->_high ? _low->_high->_action : _low->_action;
        if(action == laction && _low->_var == _var)
        {
            if(_parent->_low.get() == this)
                _parent->_low = _low;
            else
                _parent->_high = _low;
        }
    }
    if(_high && _parent && !cost_consistent && (_low == nullptr || (!_high->is_leaf() && _low->is_leaf())))
    {
        auto action = _low ? _low->_action : _action;
        auto haction = _high->_low ? _high->_low->_action : _high->_action;
        if(action == haction && _high->_var == _var)
        {
            if(_parent->_low.get() == this)
                _parent->_low = _high;
            else
                _parent->_high = _high;
        }
    }
}


bool SimpleTree::node_t::is_leaf() const {
    return _low == nullptr && _high == nullptr;
}

bool SimpleTree::node_t::can_merge(bool cost_consistent) const {
    if(_low && _low->is_leaf())
    {
        if(_low->_action != _action || (_low->_cost != _cost/* && cost_consistent*/))
            return false;
    }
    if(_high && _high->is_leaf())
    {
        if(_high->_action != _action || (_high->_cost != _cost/* && cost_consistent*/))
            return false;
    }
    assert(_low == nullptr || _high->_action == _action);
    assert(_high == nullptr || _high->_action == _action);
    return true;
}



void SimpleTree::node_t::rec_insert(double value, size_t action, std::vector<std::pair<bool,bool>>& handled, std::vector<std::pair<double, double> >& bounds, bool minimize) {
    auto nid = _var;
    assert(nid < bounds.size());
    bool reset = false;
    bool best = false;
    if( std::isnan(_cost)            ||
        (minimize && value < _cost)  ||
        (!minimize && value > _cost))
    {
        _action = action;
        _cost = value;
        best = true;
    }
    
    if(bounds[_var].second <= _limit)
    {
        if(bounds[_var].second == _limit)
            reset = handled[_var].second = true;
        if(handled[_var].first && handled[_var].second)
        {
            if(_var == handled.size()-1) {
                if(best)
                {
                    _low = nullptr;
                    _high = nullptr;
                }
                return;
            }
            ++nid;
        }
        if(is_leaf() && !best)
            return;
        if(_low == nullptr)
        {
            _low = std::make_shared<node_t>();
            _low->_var = nid;
            _low->_parent = this;
            _low->_limit = (!handled[nid].first ? bounds[nid].first : bounds[nid].second);            
        }
        _low->rec_insert(value, action, handled, bounds, minimize);
        if(reset) handled[_var].second = false;
    }
    else if(bounds[_var].first >= _limit)
    {
        if(bounds[_var].first == _limit)
            reset = handled[_var].first = true;
        if(handled[_var].first && handled[_var].second)
        {
            if(_var == handled.size()-1) {
                if(best)
                {
                    _low = nullptr;
                    _high = nullptr;
                }
                return;
            }
            ++nid;
            assert(nid < bounds.size());
        }
        if(is_leaf() && !best)
            return;
        if(_high == nullptr)
        {
            _high = std::make_shared<node_t>();
            _high->_parent = this;
            _high->_var = nid;
            _high->_limit = (!handled[nid].first ? bounds[nid].first : bounds[nid].second);
        }
        _high->rec_insert(value, action, handled, bounds, minimize);
        if(reset) handled[_var].first = false;
    }
    else
    {
        if(best && handled[_var].first && handled[_var].second && _var == handled.size()-1) {
            _low = nullptr;
            _high = nullptr;
            return;
        }
        if(is_leaf() && !best)
            return;
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

        _low->rec_insert(value, action, handled, bounds, minimize);
        _high->rec_insert(value, action, handled, bounds, minimize);
    }
}

std::ostream& SimpleTree::print(std::ostream& stream)
{
    auto dest = std::make_unique<double[]>(_statevars.size());
    stream << "{\"regressors\":\n";
    bool fe = true;
    for(size_t i = 0; i < _regressors.size(); ++i)
    {
        auto size = _regressors.unpack(i, (unsigned char*)dest.get());
        if(!fe)
            stream << ",";
        stream << "\"(";
        bool fd = true;
        for(size_t d = 0; d < size/sizeof(double); ++d)
        {
            if(!fd)
                stream << ",";
            stream << dest[d];
            fd = false;
        }
        stream << ")\":\n";
        _regressors.get_data(i).print(stream, this, 1);
        stream << "\n";
        fe = false;
    }
    stream << "}";
    return stream;
}

std::ostream& SimpleTree::node_t::print(std::ostream& out, SimpleTree* parent, size_t tabs) const {
    for(size_t i = 0; i < tabs; ++i) out << "\t";
    out << "{\"var\":" << _var << ",\"bound\":" << _limit;
    if(!std::isnan(_cost))
    {
        out << ",\"value\":" << _cost << ",\"action\":" << _action;
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

std::ostream& SimpleTree::print_c(std::ostream& stream, std::string name, size_t split) 
{
    stream << "int " << name << "(const bool* active, const double* disc, const double* vars)\n{\n";
    auto dest = std::make_unique<double[]>(_statevars.size());
    bool first = true;
    for(size_t i = 0; i < _regressors.size(); ++i)
    {
        auto size = _regressors.unpack(i, (unsigned char*)dest.get());
        bool fd = true;
        stream << "\t";
        if(!first)
            stream << " else ";
        stream << "if(";
        for(size_t d = 0; d < size/sizeof(double); ++d)
        {
            if(!fd)
                stream << " && ";
            if(d < split)
                stream << "active[" << d << "]";
            else 
                stream << "disc[" << (d-split) << "]";
            stream << " == " << dest[d];
            fd = false;
        }
        stream << ") {";
        first = false;
        _regressors.get_data(i).print_c(stream, 2);
        stream << "\t}\n";
    }
    stream << "\t else { return -1; }\n";
    stream << "}\n";
    return stream;
}

std::ostream& SimpleTree::node_t::print_c(std::ostream& stream, size_t tabs) {
    if(is_leaf())
    {
        for(size_t i = 0; i < tabs; ++i) stream << "\t";
        stream << "return " << _action << ";" << std::endl;
        return stream;
    }
    if(std::isinf(_limit))
    {
        if(_low)
            return _low->print_c(stream, tabs);
        else if(_high)
            return _high->print_c(stream, tabs);
        return stream;
    }
    for(size_t i = 0; i < tabs; ++i) stream << "\t";
    stream << "if(vars[" << _var << "] <= " << _limit << ") {\n";
    if(_low)
        _low->print_c(stream, tabs + 1);
    else
    {
        for(size_t i = 0; i < tabs + 1; ++i) stream << "\t";
        stream << "return " << _action << ";\n";
    }
    for(size_t i = 0; i < tabs; ++i) stream << "\t";
    stream << "} else {\n";
    if(_high)
        _high->print_c(stream, tabs + 1);
    else
    {
        for(size_t i = 0; i < tabs + 1; ++i) stream << "\t";
        stream << "return " << _action << ";\n";        
    }
    stream << "\n";
    for(size_t i = 0; i < tabs; ++i) stream << "\t";
    stream << "}\n";
    return stream;    
}
