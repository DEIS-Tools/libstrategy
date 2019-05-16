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
        std::cerr << "HANDLING " << it.key() << std::endl;
        for(auto iit = reg["regressor"].begin(); iit != reg["regressor"].end(); ++iit)
        {
            std::cerr << "\tHANDLING " << iit.key() << std::endl; 
            size_t action = atoi(iit.key().c_str());
            if(action >= tree._actions.size())
                throw base_error("Action-index is out of bounds");
            
            if(tree._root == nullptr)
            {
                tree._root = std::make_shared<node_t>();
                tree._root->_limit = -std::numeric_limits<double>::infinity();
                tree._root->_cost = std::numeric_limits<double>::quiet_NaN();
            }
            auto obj = iit.value();
            tree._root->merge(key, obj, action, tree._pointvars.size(), is_minimize, tree);
            tree._root = tree._root->simplify(true, tree._nodemap);
        }
        tree._root = tree._root->simplify(false, tree._nodemap);
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

void SimpleTree::node_t::merge(std::vector<double>& key, json& tree, size_t action, size_t vars, bool minimize, SimpleTree& parent) {
    
    std::vector<std::pair<double,double>> bounds;
    for(auto k : key)
        bounds.emplace_back(k,k);
    bounds.resize(vars+key.size(), std::pair<double,double>(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()));
    find_and_insert(tree, action, bounds, minimize, key.size(), parent);
}

void SimpleTree::node_t::find_and_insert(json& tree, size_t action, std::vector<std::pair<double, double> >& bounds, bool minimize, size_t keysize, SimpleTree& parent) {
    if(tree.is_object())
    {
        size_t var = tree["var"].get<size_t>() + keysize;
        double ov = tree["bound"].get<double>();
        std::swap(ov, bounds[var].second);
        find_and_insert(tree["low"], action, bounds, minimize, keysize, parent);
        std::swap(ov, bounds[var].second);
        std::swap(ov, bounds[var].first);
        find_and_insert(tree["high"], action, bounds, minimize, keysize, parent);
        std::swap(ov, bounds[var].first);
    }
    else if(tree.is_number())
    {
        // do the merge
        std::cerr << "\t\tINSERT ";
        for(auto d : bounds)
        {
            std::cerr << "[" << d.first << ", " << d.second << "]";
        }
        std::cerr << std::endl;
        insert(tree.get<double>(), action, bounds, minimize);
    }
    else
    {
        std::cerr << tree << std::endl;
        throw base_error("Leafs of trees are expected to be doubles");
    }
}

void SimpleTree::node_t::insert(double value, size_t action, std::vector<std::pair<double, double> >& bounds, bool minimize) 
{
    std::vector<std::pair<bool,bool>> handled(bounds.size());
    rec_insert(value, action, handled, bounds, minimize);
}

void SimpleTree::simplify() {
    _nodemap = ptrie::map<std::shared_ptr<node_t>>();
    _root->simplify(false, _nodemap);
}


std::shared_ptr<SimpleTree::node_t> SimpleTree::node_t::simplify(bool cost_consistent, ptrie::map<std::shared_ptr<node_t>>& nodemap) {
    if(_low) _low = _low->simplify(cost_consistent, nodemap);
    if(_high) _high = _high->simplify(cost_consistent, nodemap);
    if(std::isinf(_limit) && (_low != nullptr || _high != nullptr))
    {
        return _low == nullptr ? _high : _low;
    }
    assert(_low.get() != this);
    assert(_high.get() != this);
    if(_low && _low->is_leaf() && (_low->_cost == _cost || !cost_consistent) && _low->_action == _action) 
        _low = nullptr;
    if(_high && _high->is_leaf() && (_high->_cost == _cost || !cost_consistent) && _high->_action == _action) 
        _high = nullptr;
    if(!cost_consistent)
    {
        if( (_low == nullptr || _low->is_leaf()) &&
            (_high == nullptr || _high->is_leaf()))
        {
            auto lact = _low ? _low->_action : _action;
            auto hact = _high ? _high->_action : _action;
            if(lact == hact)
            {
                _low = nullptr;
                _high = nullptr;
                _action = lact;
                _limit = std::numeric_limits<double>::infinity();
            }
        }
    }
    if(_low && _parent && !cost_consistent && (_high == nullptr || (!_low->is_leaf() && _high->is_leaf())))
    {
        auto action = _high ? _high->_action : _action;
        auto laction = _low->_high ? _low->_high->_action : _low->_action;
        if(action == laction && _low->_var == _var)
            return _low;
    }
    if(_high && _parent && !cost_consistent && (_low == nullptr || (!_high->is_leaf() && _low->is_leaf())))
    {
        auto action = _low ? _low->_action : _action;
        auto haction = _high->_low ? _high->_low->_action : _high->_action;
        if(action == haction && _high->_var == _var)
            return _high;
    }
    if(!cost_consistent)
    {
        if(_high != nullptr && _high == _low)
            return _high;
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
        ((uint32_t*)next)[0] = _action;
        res.second = sizeof(uint32_t);
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
    bool reset_low = false;
    bool reset_high = false;
    bool best = false;
    if( std::isnan(_cost)            ||
        (minimize && value < _cost)  ||
        (!minimize && value > _cost))
    {
        _action = action;
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
        _low->rec_insert(value, action, handled, bounds, minimize);
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
        _high->rec_insert(value, action, handled, bounds, minimize);
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

        _low->rec_insert(value, action, handled, bounds, minimize);
        _high->rec_insert(value, action, handled, bounds, minimize);
    }
    if(reset_low) handled[_var].first = false;
    if(reset_high) handled[_var].second = false;
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
    stream << "\t// Depth = " << _root->depth() << std::endl;
    _root->print_c(stream, _statevars.size(), 1);
    for(size_t i = 0; i < _nodemap.size(); ++i)
    {
        auto& node = _nodemap.get_data(i);
        stream << "\tl" << node.get() << ":\n";
        node->print_c(stream, _statevars.size(), 2);
    }
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


std::ostream& SimpleTree::node_t::print_c(std::ostream& stream, size_t disc, size_t tabs) {

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
    auto type = disc > _var ? "disc" : "vars";
    auto vid = disc > _var ? _var : _var - disc;
    stream << "if(" << type << "[" << vid << "] <= " << _limit << ") ";
    if(_low)
    {
        if(_low->is_leaf() && std::isinf(_low->_limit))
            stream << "return " << _low->_action;
        else
            stream << "goto l" << _low.get();
    }
    else
        stream << "return " << _action;
    stream << ";\n";
    for(size_t i = 0; i < tabs; ++i) stream << "\t";
    stream << "else ";
    if(_high)
    {
        if(_high->is_leaf() && std::isinf(_high->_limit))
            stream << "return " << _high->_action;
        else
            stream << "goto l" << _high.get();
    }
    else
        stream << "return " << _action;
    stream << ";\n";   
    return stream;    
}
