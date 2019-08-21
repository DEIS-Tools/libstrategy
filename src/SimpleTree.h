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
    static SimpleTree parse(std::istream&, bool simplify = true, bool subsumption = true, double accuracy = 0);
    static SimpleTree parse(std::istream&, bool simplify, bool subsumption, double accuracy, std::vector<double>& exactness);
    std::ostream& print(std::ostream& stream) const;
    std::ostream& print_c(std::ostream& stream, std::string name) const;
    double value(const double* disc, const double* cont, uint32_t action) const;
    bool is_minimization() const { return _is_minimization; }
    const std::vector<std::string> &actions() const { return _actions; }
    const std::vector<std::string> &discrete_features() const { return _statevars; }
    const std::vector<std::string> &continous_features() const { return _pointvars; }
private:
    using json = nlohmann::json;
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
        std::shared_ptr<node_t> simplify(bool make_dd, ptrie::map<std::shared_ptr<node_t>>& nodemap, SimpleTree& parent);
        void subsumption_reduction(bool minimization, SimpleTree& parent);
        void action_nodes(std::vector<std::shared_ptr<node_t>>& nodes, uint32_t low, uint32_t high, uint32_t varid);
        std::pair<double,double> compute_min_max();
        bool check_tiles(node_t* start, std::vector<std::shared_ptr<node_t>>& , std::vector<std::pair<double,double>>& bounds, double val, double minval, double maxval, bool minimization, size_t offset);
        bool subsumes(std::vector<std::pair<double,double>>& bounds, std::vector<std::pair<double,double>>& obounds, double val, bool minimization, size_t offset, double& best, std::pair<double,double>& closest);
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
        std::pair<std::unique_ptr<unsigned char[]>, size_t> signature() const;
    };
        
    std::vector<std::string> _actions;
    std::vector<std::string> _statevars;
    std::vector<std::string> _pointvars;
    std::shared_ptr<node_t> _root;
    bool _is_minimization = true;
    
};

#endif /* SIMPLETREE_H */

