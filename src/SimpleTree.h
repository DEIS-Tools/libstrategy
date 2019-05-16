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

#include <ptrie_map.h>
#include <nlohmann/json.hpp>

class SimpleTree {
public:
    SimpleTree(const SimpleTree& orig) = default;
    virtual ~SimpleTree() = default;
    static SimpleTree parse(std::istream&);
    std::ostream& print(std::ostream& stream);
    std::ostream& print_c(std::ostream& stream, std::string name, size_t split);
    void simplify();
private:
    using json = nlohmann::json;
    SimpleTree() = default;
    
    static std::vector<double> parse_key(const std::string& key);
    
    struct node_t : public std::enable_shared_from_this<node_t> {
        uint32_t _var = 0;
        uint32_t _action = 0;
        double _limit = -std::numeric_limits<double>::infinity();
        double _cost = std::numeric_limits<double>::quiet_NaN();
        std::shared_ptr<node_t> _low;
        std::shared_ptr<node_t> _high;
        node_t* _parent;
        void merge(std::vector<double>& key, json& tree, size_t action, size_t vars, bool minimize, SimpleTree& parent);
        void find_and_insert(json& tree, size_t action, std::vector<std::pair<double,double>>& bounds, bool minimize, size_t keysize, SimpleTree& parent);
        void insert(double value, size_t action, std::vector<std::pair<double,double>>& bounds, bool minimize);
        void rec_insert(double value, size_t action, std::vector<std::pair<bool,bool>>& handled, std::vector<std::pair<double, double> >& bounds, bool minimize);
        std::ostream& print(std::ostream& out, SimpleTree* parent, size_t tabs = 0) const;
        bool is_leaf() const;
        bool can_merge(bool cost_consistent) const;
        std::shared_ptr<node_t> simplify(bool cost_consistent, ptrie::map<std::shared_ptr<node_t>>& nodemap);
        std::ostream& print_c(std::ostream& stream, size_t disc, size_t tabs = 0);
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

