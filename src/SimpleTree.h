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
    
    struct node_t {
        uint32_t _var = 0;
        uint32_t _action = 0;
        double _limit = -std::numeric_limits<double>::infinity();
        double _cost = std::numeric_limits<double>::quiet_NaN();
        std::shared_ptr<node_t> _low;
        std::shared_ptr<node_t> _high;
        node_t* _parent;
        void merge(json tree, size_t action, size_t vars, bool minimize);
        void find_and_insert(json tree, size_t action, std::vector<std::pair<double,double>>& bounds, bool minimize);
        void insert(double value, size_t action, std::vector<std::pair<double,double>>& bounds, bool minimize);
        void rec_insert(double value, size_t action, std::vector<std::pair<bool,bool>>& handled, std::vector<std::pair<double, double> >& bounds, bool minimize);
        std::ostream& print(std::ostream& out, SimpleTree* parent, size_t tabs = 0) const;
        bool is_leaf() const;
        bool can_merge(bool cost_consistent) const;
        void simplify(bool cost_consistent);
        std::ostream& print_c(std::ostream& stream, size_t tabs = 0);
    };
        
    std::vector<std::string> _actions;
    std::vector<std::string> _statevars;
    std::vector<std::string> _pointvars;
    ptrie::map<node_t> _regressors;
    
    
};

#endif /* SIMPLETREE_H */

