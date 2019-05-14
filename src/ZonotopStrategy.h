/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ZonotopStrategy.h
 * Author: petko
 *
 * Created on December 13, 2018, 2:54 PM
 */

#ifndef ZONOTOPSTRATEGY_H
#define ZONOTOPSTRATEGY_H

#include <ptrie_stable.h>
#include <istream>
#include <memory>

#include "SimpleTree.h"
#include "errors.h"

class ZonotopStrategy {
public:
    ZonotopStrategy(const ZonotopStrategy& orig) = default;
    virtual ~ZonotopStrategy() = default;
    static ZonotopStrategy parse(std::istream&);
    int num_patterns() const;
    int max_pattern_length() const;
    void active(const double* sample, bool* write) const;
    int get_pattern(int el, int* write);
    double get_min(size_t dimen) const;
    double get_max(size_t dimen) const;
    std::ostream& print(std::ostream&);
    std::ostream& print_c(std::ostream& stream, std::string function_name);
private:
    ZonotopStrategy() = default;
    ptrie::set_stable<> _states;
    size_t _max_length = 0;
    struct node_t {
        uint32_t _varid = 0;        
        double _limit = 0;
        std::shared_ptr<node_t> _low = nullptr;
        std::shared_ptr<node_t> _high = nullptr;
        std::vector<size_t> _low_patterns;
        std::vector<size_t> _high_patterns;
        node_t* _parent;
        std::ostream& print(std::ostream& out, ZonotopStrategy* parent, size_t tabs = 0) const;
        std::ostream& print_c(std::ostream& out, size_t tabs = 0) const;
        double get_min(size_t dimen) const;
        double get_max(size_t dimen) const;
    };
    std::shared_ptr<node_t> _root;
    struct bound_t {
        double _lower = 0;
        double _upper = 0;
        bound_t() = default;
        bound_t(const bound_t&) = default;
        bound_t(bound_t&&) = default;
        bound_t(double l, double u) 
        : _lower(l), _upper(u) {}
    };
    static std::vector<bound_t> get_bounds(const std::string& zonotop, size_t& vars);
    bool add(std::vector<bound_t>& bounds, size_t state);
    bool rec_insert(node_t* node, std::vector<bound_t>& bounds, std::vector<std::pair<bool,bool>>& handled, size_t state);
    void try_merge(node_t* node, size_t state);
};


#endif /* ZONOTOPSTRATEGY_H */

