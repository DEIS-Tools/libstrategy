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
 * File:   ZonotopStrategy.h
 * Author: Peter G. Jensen
 *
 * Created on December 13, 2018, 2:54 PM
 */

#ifndef ZONOTOPSTRATEGY_H
#define ZONOTOPSTRATEGY_H

#include <ptrie/ptrie_stable.h>
#include <istream>
#include <memory>

#include "SimpleTree.h"
#include "errors.h"

class ZonotopStrategy {
public:
    ZonotopStrategy(ZonotopStrategy&&) = default;
    virtual ~ZonotopStrategy() = default;
    ZonotopStrategy& operator=(ZonotopStrategy&&) = default;
    static ZonotopStrategy parse(std::istream&);
    int num_patterns() const;
    int max_pattern_length() const;
    void active(const double* sample, bool* write) const;
    int get_pattern(int el, int* write);
    double get_min(size_t dimen) const;
    double get_max(size_t dimen) const;
    std::ostream& print(std::ostream&) const;
    std::ostream& print_c(std::ostream& stream, std::string function_name) const;
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
        std::ostream& print(std::ostream& out, const ZonotopStrategy* parent, size_t tabs = 0) const;
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

