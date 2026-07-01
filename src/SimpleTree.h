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

#include <ptrie/ptrie_map.h>  // includes iostream :-(

#include <iosfwd>
#include <memory>
#include <vector>

class SimpleTree
{
public:
    SimpleTree(const SimpleTree&) = delete;
    SimpleTree(SimpleTree&&) noexcept = default;
    SimpleTree& operator=(const SimpleTree&) = delete;
    SimpleTree& operator=(SimpleTree&&) noexcept;
    ~SimpleTree();
    static SimpleTree parse(std::istream&, bool simplify = false, bool subsumption = false, double accuracy = 0);
    static SimpleTree parse(std::istream&, bool simplify, bool subsumption, double accuracy,
                            std::vector<double>& exactness);

    std::ostream& print(std::ostream& os) const;
    std::ostream& print_c(std::ostream& os, const std::string& name) const;
    double value(const double* disc, const double* cont, uint32_t action) const;
    bool is_minimization() const { return _is_minimization; }
    const std::vector<std::string>& actions() const;
    const std::vector<std::string>& discrete_features() const { return _statevars; }
    const std::vector<std::string>& continous_features() const { return _pointvars; }

private:
    struct node_t;
    using node_ptr = std::shared_ptr<node_t>;
    struct signature_t
    {
        uint32_t _var{0};
        double _limit{0};
        node_t* _low{nullptr};
        node_t* _high{nullptr};
    } __attribute__((packed));
    friend struct ptrie::byte_iterator<signature_t>;

    using nodemap_t = ptrie::map<signature_t, node_ptr>;
    SimpleTree() = default;

    std::vector<std::string> _actions;
    std::vector<std::string> _statevars;
    std::vector<std::string> _pointvars;
    node_ptr _root;
    bool _is_minimization = true;
};

template <>
struct ptrie::byte_iterator<SimpleTree::signature_t>
{
    static uchar& access(SimpleTree::signature_t* data, size_t id) { return reinterpret_cast<uchar*>(data)[id]; }

    static const uchar& const_access(const SimpleTree::signature_t* data, size_t id)
    {
        return reinterpret_cast<const uchar*>(data)[id];
    }

    static constexpr size_t element_size()
    {
        constexpr auto member_size = sizeof(SimpleTree::signature_t::_var) + sizeof(SimpleTree::signature_t::_limit) +
                                     sizeof(SimpleTree::signature_t::_high) + sizeof(SimpleTree::signature_t::_low);
        static_assert(sizeof(SimpleTree::signature_t) == member_size, "tightly packed struct");
        return sizeof(SimpleTree::signature_t);
    }

    [[deprecated]] static constexpr bool continious() { return true; }
    static constexpr bool continuous() { return true; }

    // add read_blob, write_blob
};
#endif /* SIMPLETREE_H */
