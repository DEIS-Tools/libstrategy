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
 * File:   SimpleTree.h
 * Author: Peter G. Jensen <root@petergjoel.dk>
 *
 * Created on May 9, 2019, 10:21 PM
 */

#ifndef SIMPLETREE_H
#define SIMPLETREE_H

#include <istream>
#include <memory>
#include <vector>
#include <unordered_set>
#include <set>

#include <ptrie/ptrie_map.h>
#include <nlohmann/json.hpp>

class SimpleTree {
public:
    SimpleTree(const SimpleTree& orig) = default;
    virtual ~SimpleTree() = default;
    static SimpleTree parse(std::istream&, bool simplify = false, bool subsumption = false, double accuracy = 0);
    static SimpleTree parse(std::istream&, bool simplify, bool subsumption, double accuracy, std::vector<double>& exactness);
    std::ostream& print(std::ostream& stream) const;
    std::ostream& print_c(std::ostream& stream, std::string name) const;
    double value(const double* disc, const double* cont, uint32_t action) const;
    bool is_minimization() const { return _is_minimization; }
    const std::vector<std::string> &actions() const;
    const std::vector<std::string> &discrete_features() const { return _statevars; }
    const std::vector<std::string> &continous_features() const { return _pointvars; }
private:
    
    using json = nlohmann::json;
    struct node_t;
    struct signature_t {
        uint32_t _var;
        double _limit;
        node_t* _low;
        node_t* _high;
        signature_t(const SimpleTree::node_t&);
    } __attribute__((packed));
    friend struct ptrie::byte_iterator<signature_t>;

    using nodemap_t = ptrie::map<signature_t,std::shared_ptr<node_t>>;
    SimpleTree() = default;
    
    static std::vector<double> parse_key(const std::string& key);
    
    struct node_t : public std::enable_shared_from_this<node_t> {
        uint32_t _var = std::numeric_limits<uint32_t>::max();
        double _limit = -std::numeric_limits<double>::infinity();
        double _cost = std::numeric_limits<double>::infinity();
        std::pair<double,double> _cost_bounds;
        std::shared_ptr<node_t> _low;
        std::shared_ptr<node_t> _high;
        node_t* _parent;
        void insert(std::vector<double>& key, json& tree, size_t action, SimpleTree& parent, size_t prefix, bool minimize, double accuracy, std::vector<double>& exactness);
        std::ostream& print(std::ostream& out, size_t tabs = 0) const;
        bool is_leaf() const;
        std::shared_ptr<node_t> simplify(bool make_dd, nodemap_t& nodemap, SimpleTree& parent);
        void subsumption_reduction(bool minimization, SimpleTree& parent);
        void action_nodes(std::vector<std::shared_ptr<node_t>>& nodes, uint32_t low, uint32_t high, uint32_t varid);
        std::pair<double,double> compute_min_max();
        bool check_tiles(node_t* start, std::vector<std::shared_ptr<node_t>>& , std::vector<std::pair<double,double>>& bounds, double val, double minval, double maxval, bool minimization, size_t offset);
        bool subsumes(const std::vector<std::pair<double,double>>& bounds, std::vector<std::pair<double,double>>& obounds, const double val, const bool minimization, size_t offset, double& best, std::pair<double,double>& closest);
        void get_ranks(std::set<std::pair<double, node_t*>>& values, node_t* start);
        void set_ranks(std::unordered_map<double,double>& values);
        std::ostream& print_c(std::ostream& stream, size_t disc, std::unordered_set<const node_t*>& printed, size_t tabs = 0) const;
        std::ostream& print_c_nested(std::ostream& stream, size_t disc, size_t tabs, std::vector<const node_t*>& toprint, const std::shared_ptr<node_t>& node) const;
        size_t depth() const;
        double value(const double* disc, const double* cont, uint32_t action, size_t ndisc) const;
        void consistent(size_t) const;
        bool cost_intersect(const node_t& other) const;
        double midcost(const node_t& other, double minval, double maxval) const;
        bool operator==(const node_t& other) const
        {
            if(_limit != other._limit)
                return false;
            if(_var != other._var)
                return false;
            if(_low != other._low)
                return false;
            if(_high != other._high)
                return false;
            return true;
        }
        bool operator<(const node_t& other) const
        {
            if(_limit != other._limit)
                return _limit < other._limit;
            if(_var != other._var)
                return _var < other._var;
            if(_low != other._low)
                return _low < other._low;
            return _high < other._high;
        }
        std::shared_ptr<node_t>& operator[](bool b) {
            return b ? _high : _low;
        }
    };
        
    std::vector<std::string> _actions;
    std::vector<std::string> _statevars;
    std::vector<std::string> _pointvars;
    std::shared_ptr<node_t> _root;
    bool _is_minimization = true;
};

namespace ptrie{
    template<>
    struct byte_iterator<SimpleTree::signature_t> {
        static uchar& access(SimpleTree::signature_t* data, size_t id)
        {
            return ((uchar*)data)[id];
        }
        
        static const uchar& const_access(const SimpleTree::signature_t* data, size_t id)
        {
            return ((const uchar*)data)[id];
        }
        
        static constexpr size_t element_size()
        {
            return sizeof(size_t)*2+sizeof(double)+sizeof(uint32_t);
        }
        
        static constexpr bool continious()
        {
            return true;
        }
        
        // add read_blob, write_blob
    };
}
#endif /* SIMPLETREE_H */

