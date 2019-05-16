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

#include <ptrie_map.h>
#include <nlohmann/json.hpp>

class SimpleTree {
public:
    SimpleTree(const SimpleTree& orig) = default;
    virtual ~SimpleTree() = default;
    static SimpleTree parse(std::istream&);
    std::ostream& print(std::ostream& stream);
    std::ostream& print_c(std::ostream& stream, std::string name, size_t split);
private:
    using json = nlohmann::json;
    SimpleTree() = default;
    
    static std::vector<double> parse_key(const std::string& key);
    
    struct node_t : public std::enable_shared_from_this<node_t> {
        uint32_t _var = std::numeric_limits<uint32_t>::max();
        double _limit = -std::numeric_limits<double>::infinity();
        double _cost = std::numeric_limits<double>::infinity();
        double _mincost = std::numeric_limits<double>::infinity();
        double _maxcost = -std::numeric_limits<double>::infinity();
        std::shared_ptr<node_t> _low;
        std::shared_ptr<node_t> _high;
        node_t* _parent;
        void insert(std::vector<double>& key, json& tree, size_t action, SimpleTree& parent, size_t prefix, bool minimize);
        std::ostream& print(std::ostream& out, size_t tabs = 0) const;
        bool is_leaf() const;
        std::shared_ptr<node_t> simplify(bool make_dd, ptrie::map<std::shared_ptr<node_t>>& nodemap);
        void subsumption_reduction(bool minimization, SimpleTree& parent);
        void action_nodes(std::vector<std::shared_ptr<node_t>>& nodes, uint32_t low, uint32_t high, uint32_t varid);
        std::pair<double,double> compute_min_max();
        void check_tiles(node_t* start, std::vector<std::shared_ptr<node_t>>& , std::vector<std::pair<double,double>>& bounds, double val, bool minimization, size_t offset);
        bool subsumes(std::vector<std::pair<double,double>>& bounds, double val, bool minimization, size_t offset);
        void get_ranks(std::set<std::pair<double, node_t*>>& values, node_t* start);
        void set_ranks(std::unordered_map<double,double>& values);
        std::ostream& print_c(std::ostream& stream, size_t disc, std::unordered_set<node_t*>& printed, size_t tabs = 0);
        size_t depth() const;        
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
    std::shared_ptr<node_t> _root = nullptr;
    ptrie::map<std::shared_ptr<node_t>> _nodemap;
    
};

#endif /* SIMPLETREE_H */

