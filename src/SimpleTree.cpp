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
 * File:   SimpleTree.cpp
 * Author: Peter G. Jensen <root@petergjoel.dk>
 *
 * Created on May 9, 2019, 10:21 PM
 */

#include "SimpleTree.h"
#include "errors.h"
#include "utilities.hpp"

#include <nlohmann/json.hpp>

#include <iostream>
#include <unordered_set>
#include <vector>

using json = nlohmann::json;

const std::vector<std::string> &SimpleTree::actions() const { return _actions; }

SimpleTree SimpleTree::parse(std::istream &input, bool simplify,
                             bool subsumption, double accuracy) {
  auto empty = std::vector<double>{};
  return parse(input, simplify, subsumption, accuracy, empty);
}

SimpleTree SimpleTree::parse(std::istream &input,
                             bool simplify [[maybe_unused]], bool subsumption,
                             double accuracy, std::vector<double> &exactness) {
  auto raw = json::parse(input);
  if (!raw.is_object())
    throw base_error("Input JSON not well formatted");

  if (!raw["version"].is_number_float() || raw["version"].get<double>() != 1.0)
    throw base_error("Input version not supported");

  if (!raw["type"].is_string() ||
      raw["type"].get<std::string>() != "state->regressor")
    throw base_error("Input type not supported");

  if (!raw["representation"].is_string() ||
      raw["representation"].get<std::string>() != "map")
    throw base_error("Input representation not supported");

  if (!raw["actions"].is_object())
    throw base_error("Expected actions object");
  // store actions
  SimpleTree tree;
  if (tree._root == nullptr) {
    tree._root = std::make_shared<node_t>();
    tree._root->_limit = -std::numeric_limits<double>::infinity();
    tree._root->_cost = std::numeric_limits<double>::infinity();
    tree._root->_var = std::numeric_limits<uint32_t>::max();
  }
  for (const auto &[act_key, act_value] : raw["actions"].items()) {
    auto id = std::stoi(act_key);
    if (id >= static_cast<int>(tree._actions.size()))
      tree._actions.resize(id + 1);
    tree._actions[id] = act_value.get<std::string>();
  }

  if (!raw["statevars"].is_array())
    throw base_error("Expected statevars as array");
  tree._statevars = raw["statevars"].get<std::vector<std::string>>();

  if (!raw["pointvars"].is_array())
    throw base_error("Expected pointvars as array");
  tree._pointvars = raw["pointvars"].get<std::vector<std::string>>();

  if (!raw["regressors"].is_object())
    throw base_error("Expected regressors as object");
  auto minim = -1;
  auto inf = json("inf"); // json string (whereas {"inf"} is json array!)
  auto first_element = true;
  for (const auto &[raw_key, raw_value] : raw["regressors"].items()) {
    auto key = parse_key(raw_key);
    if (key.size() != tree._statevars.size() &&
        (!tree._statevars.empty() || key[0] != 1.0)) {
      std::cerr << raw_key << std::endl;
      std::cerr << key.size() << std::endl;
      std::cerr << tree._statevars.size() << std::endl;
      throw base_error("Expected cardinality of key does not match cardinality "
                       "of statevars");
    }
    if (tree._statevars.empty())
      key.clear();

    if (!raw_value.is_object())
      throw base_error("Expected regressor-element as object");
    if (raw_value["type"].get<std::string>() != "act->point->val")
      throw base_error("Expected regressor-element of type act->point->val");
    if (raw_value["representation"].get<std::string>() != "simpletree")
      throw base_error("Expected regressor-element represented as simpletree");
    bool is_minimize = raw_value["minimize"].is_boolean()
                           ? raw_value["minimize"].get<bool>()
                           : raw_value["minimize"].get<int32_t>();
    if (!raw_value["regressor"].is_object())
      throw base_error("Expected regressor-sub-element as an object");
    if (minim == -1)
      minim = is_minimize;
    if (minim != is_minimize)
      throw base_error(
          "Expected all sub-regressors to have same minimization flag");
    if (first_element)
      tree._root->_cost =
          (is_minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
    first_element = false;
    // make sure all actions are mapped initially
    for (size_t i = 0; i < tree._actions.size(); ++i)
      tree._root->insert(key, inf, i, tree, 0, is_minimize, 0, exactness);
    for (const auto &[reg_key, reg_value] : raw_value["regressor"].items()) {
      auto action = std::stoi(reg_key);
      if (action >= static_cast<int>(tree._actions.size()))
        throw base_error("Action-index is out of bounds");

      if (tree._root == nullptr) {
        tree._root = std::make_shared<node_t>();
        tree._root->_limit = -std::numeric_limits<double>::infinity();
        tree._root->_cost = std::numeric_limits<double>::infinity();
        tree._root->_cost =
            (is_minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
        tree._root->_var = std::numeric_limits<uint32_t>::max();
      }

      tree._root->insert(key, reg_value, action, tree, 0, is_minimize, accuracy,
                         exactness);
    }
  }
  tree._is_minimization = minim != 0;
  if (subsumption)
    tree._root->subsumption_reduction(minim, tree);
  /*if(simplify) // disabled for now
  {
      nodemap_t nodemap;
      if(tree._root)
          tree._root = tree._root->simplify(true, nodemap, tree);
      if(tree._root)
          tree._root->_parent = nullptr;
  }*/
  return tree;
}

std::pair<double, double> SimpleTree::node_t::compute_min_max() const {
  double mincost = std::numeric_limits<double>::infinity();
  double maxcost = -mincost;
  if (_low) {
    const auto [r_min, r_max] = _low->compute_min_max();
    mincost = std::min(r_min, mincost);
    maxcost = std::max(r_max, maxcost);
  }
  if (_high) {
    const auto [r_min, r_max] = _high->compute_min_max();
    mincost = std::min(r_min, mincost);
    maxcost = std::max(r_max, maxcost);
  }

  if (is_leaf() && !std::isinf(_cost)) {
    mincost = std::min(_cost, mincost);
    maxcost = std::max(_cost, maxcost);
  }
  return std::make_pair(mincost, maxcost);
}

void SimpleTree::node_t::action_nodes(std::vector<node_ptr> &nodes,
                                      uint32_t low, uint32_t high,
                                      uint32_t varid) {
  if (_var != varid) {
    if (!is_leaf() || !(std::isnan(_cost) || std::isinf(_cost))) {
      auto ptr = shared_from_this();
      auto lb = std::lower_bound(nodes.begin(), nodes.end(), ptr);
      if (lb == nodes.end() || lb->get() != this)
        nodes.insert(lb, ptr);
    }
  } else {
    if (_low)
      _low->action_nodes(nodes, low, _limit, varid);
    if (_high)
      _high->action_nodes(nodes, _limit + 1, high, varid);
  }
}

void SimpleTree::node_t::subsumption_reduction(bool minimization,
                                               SimpleTree &parent) {
  if (_var < parent._statevars.size()) {
    if (_low)
      _low->subsumption_reduction(minimization, parent);
    if (_high)
      _high->subsumption_reduction(minimization, parent);
  } else {
    if (_var != parent._statevars.size()) {
      return;
    }
    size_t prev_vals = std::numeric_limits<size_t>::max();
    size_t p = 0;
    while (true) {
      ++p;
      bool outer = false;
      std::vector<node_ptr> nodes;
      action_nodes(nodes, 0, parent._actions.size() - 1,
                   parent._statevars.size());
      auto val = std::numeric_limits<double>::infinity();
      if (!minimization)
        val *= -1;
      auto best = val;
      auto worst = -val;
      for (const auto &n : nodes) {
        auto [min, max] = n->compute_min_max();
        if (std::isinf(min))
          continue;
        assert(!std::isinf(max));
        if (minimization) {
          val = std::min(val, max);
          best = std::min(best, min);
          worst = std::max(worst, max);
        } else {
          val = std::max(val, min);
          best = std::max(best, max);
          worst = std::min(worst, min);
        }
      }
      std::vector<std::pair<double, double>> bounds(
          parent._pointvars.size(),
          std::make_pair(-std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity()));
      for (const auto &n : nodes) {
        outer |= n->check_tiles(n.get(), nodes, bounds, val, best, worst,
                                minimization, parent._statevars.size() + 1);
      }
      std::set<std::pair<double, node_t *>> values;
      {
        auto nodes = std::vector<node_ptr>(parent._actions.size());
        action_nodes(nodes, 0, parent._actions.size() - 1,
                     parent._statevars.size());
        for (auto &n : nodes) {
          if (n == nullptr)
            continue;
          get_ranks(values, n.get());
        }
        std::unordered_map<double, double> replace;
        for (auto &e : values) {
          if (replace.count(e.first) > 0)
            continue;
          double id = replace.size();
          replace[e.first] = id;
          //            std::cerr << e.first << " : " << e.second << std::endl;
        }
        set_ranks(replace);
      }
      if (!outer && prev_vals <= values.size())
        break;
      prev_vals = values.size();
    }
  }
}

void SimpleTree::node_t::set_ranks(std::unordered_map<double, double> &values) {
  if (is_leaf()) {
    if (std::isinf(_cost) || std::isnan(_cost))
      return;
    _cost = values[_cost];
  } else {
    if (_low)
      _low->set_ranks(values);
    if (_high)
      _high->set_ranks(values);
  }
}

void SimpleTree::node_t::get_ranks(
    std::set<std::pair<double, node_t *>> &values, node_t *start) {
  if (is_leaf()) {
    if (std::isinf(_cost) || std::isnan(_cost))
      return;
    values.emplace(_cost, start);
  } else {
    if (_low)
      _low->get_ranks(values, start);
    if (_high)
      _high->get_ranks(values, start);
  }
}

bool SimpleTree::node_t::subsumes(
    const std::vector<std::pair<double, double>> &bounds,
    std::vector<std::pair<double, double>> &obounds, const double val,
    const bool minimization, size_t offset, double &best,
    std::pair<double, double> &closest) const {
  if (is_leaf()) {
    // TODO: continue in other trees here, maybe combination is
    // enough to prove subsumption!
    if (_cost <= val)
      closest.first = std::max(_cost, closest.first);
    if (_cost >= val)
      closest.second = std::min(_cost, closest.second);
    if (minimization)
      best = std::min(_cost, best);
    else
      best = std::max(_cost, best);
    if (minimization && _cost <= val)
      return true;
    if (!minimization && _cost >= val)
      return true;
    return false;
  } else {
    /*if(minimization && _maxcost <= val)
    {
        best = std::min(_maxcost, best);
        return true;
    }
    else if(!minimization && _mincost >= val)
    {
        best = std::max(_mincost, best);
        return true;
    }*/
    assert(std::isinf(_cost));
    assert(_var >= offset);
    assert(bounds[_var - offset].first <= bounds[_var - offset].second);
    double org = _limit;
    std::swap(obounds[_var - offset].second, org);
    bool res = true;
    if (bounds[_var - offset].first <= _limit && _low)
      if (!_low->subsumes(bounds, obounds, val, minimization, offset, best,
                          closest))
        res = false;
    std::swap(obounds[_var - offset].second, org);
    std::swap(obounds[_var - offset].first, org);
    if (bounds[_var - offset].second > _limit && _high)
      if (!_high->subsumes(bounds, obounds, val, minimization, offset, best,
                           closest))
        res = false;
    std::swap(obounds[_var - offset].first, org);
    return res;
  }
}

bool SimpleTree::node_t::check_tiles(
    node_t *start, std::vector<node_ptr> &nodes,
    std::vector<std::pair<double, double>> &bounds, double val, double best_val,
    double worst_val, bool minimization, size_t offset) {
  auto obounds = bounds;
  if (is_leaf()) {
    double best =
        (minimization ? 1 : -1) * std::numeric_limits<double>::infinity();
    _cost_bounds.first = -std::numeric_limits<double>::infinity();
    _cost_bounds.second = std::numeric_limits<double>::infinity();
    if (_cost ==
        (minimization ? 1 : -1) * std::numeric_limits<double>::infinity())
      return false;
    for (auto &n : nodes) {
      if (n.get() == start || n == nullptr)
        continue;
      if (n->subsumes(bounds, obounds, _cost, minimization, offset, best,
                      _cost_bounds)) {
        _cost =
            (minimization ? 1 : -1) * std::numeric_limits<double>::infinity();
        _cost_bounds.second = std::numeric_limits<double>::infinity();
        return true;
      }
      assert(best <= best_val || minimization);  // doing max
      assert(best >= best_val || !minimization); // doing min
    }
    assert(best <= best_val || minimization);  // doing max
    assert(best >= best_val || !minimization); // doing min
    if ((minimization && best >= _cost) || (!minimization && best <= _cost)) {
      // std::cerr << "BEST " << _cost << " MV " << minval << " BND " <<
      // _cost_bounds.first << ": " << _cost_bounds.second << std::endl;
      _cost = best_val;
    } else {
      if (!std::isinf(_cost_bounds.first)) {
        _cost = std::nextafter(_cost_bounds.first,
                               std::numeric_limits<double>::infinity());
        if (!std::isinf(_cost_bounds.second) &&
            std::ceil(_cost) < _cost_bounds.second)
          _cost = std::ceil(_cost);
      } else if (!std::isinf(_cost_bounds.second)) {
        _cost = std::nextafter(_cost_bounds.second,
                               -std::numeric_limits<double>::infinity());
      }
    }
    return false;
  }

  if (!is_leaf()) {
    double org = _limit;
    auto &bnd = bounds[_var - offset];
    node_ptr switchnode = nullptr;
    if (bnd.second <= _limit) {
      switchnode = _low;
    }
    if (bnd.first > _limit) {
      switchnode = _high;
    }
    /*        std::cerr << "HANDLE " << std::endl;
            print(std::cerr);
            std::cerr << std::endl;*/
    bool res = false;
    if (switchnode) {
      _var = switchnode->_var;
      _limit = switchnode->_limit;
      _cost = switchnode->_cost;
      _low = switchnode->_low;
      _high = switchnode->_high;
      return check_tiles(start, nodes, bounds, val, best_val, worst_val,
                         minimization, offset);
    } else {
      std::swap(org, bnd.second);
      if (_low)
        res |= _low->check_tiles(start, nodes, bounds, val, best_val, worst_val,
                                 minimization, offset);
      std::swap(org, bnd.second);
      std::swap(org, bnd.first);
      if (_high)
        res |= _high->check_tiles(start, nodes, bounds, val, best_val,
                                  worst_val, minimization, offset);
      std::swap(org, bnd.first);
    }
    _cost_bounds.first =
        std::max(_low->_cost_bounds.first, _high->_cost_bounds.first);
    _cost_bounds.second =
        std::min(_low->_cost_bounds.second, _high->_cost_bounds.second);
    // if both are leafs and can be merged.
    if (_low->is_leaf() && _high->is_leaf() && _low->cost_intersect(*_high) &&
        std::isinf(_low->_cost) == std::isinf(_high->_cost)) {
      /*std::cerr << "[" << _low->_cost_bounds.first << ", " <<
      _low->_cost_bounds.second << "]" << std::endl; std::cerr << "[" <<
      _high->_cost_bounds.first << ", " << _high->_cost_bounds.second << "]" <<
      std::endl; std::cerr << "[" << _cost_bounds.first << ", " <<
      _cost_bounds.second << "]" << std::endl; std::cerr << "FIXING " <<
      std::endl; print(std::cerr); std::cerr << std::endl;*/
      _cost = _low->midcost(*_high, best_val, worst_val);
      // std::cerr << "NC " << _cost << std::endl;
      _low = nullptr;
      _high = nullptr;
      _var = std::numeric_limits<decltype(_var)>::max();
      _limit = std::numeric_limits<double>::infinity();
      return true;
      // std::cerr << "REMO 2" << std::endl;
    }

    {
      // check if one is leaf and sibiling has leaf on same side which
      // can be merged:
      double nc = 0;
      node_ptr node, other;
      if (!_high->is_leaf() && _high->_var == _var &&
          _high->_low->cost_intersect(*_low) && _low->is_leaf() &&
          _high->_low->is_leaf() &&
          std::isinf(_low->_cost) == std::isinf(_high->_low->_cost)) {
        node = _high;
        other = _high->_low;
        nc = other->midcost(*_low, best_val, worst_val);
        assert(!std::isnan(nc));
      } else if (!_low->is_leaf() && _low->_var == _var &&
                 _low->_high->cost_intersect(*_high) && _high->is_leaf() &&
                 _low->_high->is_leaf() &&
                 std::isinf(_high->_cost) == std::isinf(_low->_high->_cost)) {
        node = _low;
        other = _low->_high;
        nc = other->midcost(*_high, best_val, worst_val);
        assert(!std::isnan(nc));
      }

      if (node) {
        _var = node->_var;
        _limit = node->_limit;
        _low = node->_low;
        _high = node->_high;
        other->_cost = nc;
        // we could compute the values on the fly; but no time currently
        check_tiles(start, nodes, bounds, val, best_val, worst_val,
                    minimization, offset);
        return true;
      }

      if (!_low->is_leaf() && !_high->is_leaf() && _low->_var == _high->_var &&
          _low->_var == _var && _low->_high->is_leaf() &&
          _high->_low->is_leaf() && _low->_high->cost_intersect(*_high->_low) &&
          std::isinf(_low->_high->_cost) == std::isinf(_high->_low->_cost)) {
        assert(_high->_low->is_leaf());
        assert(_low->_high->is_leaf());
        _high->_low->_cost =
            _low->_high->midcost(*_high->_low, best_val, worst_val);
        _var = _low->_var;
        _limit = _low->_limit;
        _low = _low->_low;
        // we could compute the values on the fly; but no time currently
        check_tiles(start, nodes, bounds, val, best_val, worst_val,
                    minimization, offset);
        return true;
      }
    }
  }
  return false;
}

bool SimpleTree::node_t::cost_intersect(const node_t &other) const {
  return _cost_bounds.first < other._cost_bounds.second &&
         other._cost_bounds.first < _cost_bounds.second;
}

double SimpleTree::node_t::midcost(const node_t &other, double minval,
                                   double maxval) const {
  if (_cost == other._cost)
    return _cost;
  assert(std::isinf(_cost) == std::isinf(other._cost));
  auto low = std::min(other._cost_bounds.second, _cost_bounds.second);
  auto high = std::max(other._cost_bounds.first, _cost_bounds.first);
  double c = std::ceil(low);
  if (std::isinf(low))
    c = minval;
  else if (std::isinf(high))
    c = maxval;
  else {
    c = std::round(low + high / 2);
    if (c <= low || c >= high)
      c = (low + high) / 2;
  }
  return c;
}

void SimpleTree::node_t::insert(std::vector<double> &key, json &tree,
                                size_t action, SimpleTree &parent,
                                size_t prefix, bool minimize, double accuracy,
                                std::vector<double> &exactness) {
  if ((tree.is_number() || tree.is_string()) && key.size() < prefix) {
    if (is_leaf()) {
      auto nc = std::numeric_limits<double>::infinity();
      if (!minimize)
        nc *= -1;
      if (tree.is_number())
        nc = tree.get<double>();
      if (minimize)
        _cost = std::isinf(_cost) ? nc : std::min(nc, _cost);
      else
        _cost = std::isinf(_cost) ? nc : std::max(nc, _cost);
      if (accuracy != 0)
        _cost = std::round(_cost / accuracy) * accuracy;
    } else {
      assert(false);
      _low->insert(key, tree, action, parent, prefix, minimize, accuracy,
                   exactness);
      _high->insert(key, tree, action, parent, prefix, minimize, accuracy,
                    exactness);
    }
    consistent(key.size() + 1);
    return;
  } else if (is_leaf()) {
    if (std::isinf(_limit)) {
      // assert(std::isinf(_cost));
      assert(_var == std::numeric_limits<uint32_t>::max());
      _low = std::make_shared<node_t>();
      _high = std::make_shared<node_t>();
      _low->_cost = _high->_cost =
          (minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
      if (!std::isinf(_cost)) {
        std::swap(_low->_cost, _cost);
        _high->_cost = _low->_cost;
      }

      _low->_parent = _high->_parent = this;
      if (prefix == key.size()) {
        _var = prefix;
        _limit = action;
      } else if (prefix < key.size()) {
        _var = prefix;
        _limit = key[prefix];
      } else {
        const auto &tree_var = tree["var"];
        _var = key.size() + 1 + tree_var.get<uint32_t>();
        auto ivar = tree_var.get<uint32_t>();
        _limit = tree["bound"].get<double>();
        if (exactness.size() > ivar && exactness[ivar] != 0)
          _limit = std::round(_limit / exactness[ivar]) * exactness[ivar];
      }
      consistent(key.size() + 1);
      if (prefix <= key.size())
        _low->insert(key, tree, action, parent, prefix + 1, minimize, accuracy,
                     exactness);
      else {
        _low->insert(key, tree["low"], action, parent, prefix + 1, minimize,
                     accuracy, exactness);
        _high->insert(key, tree["high"], action, parent, prefix + 1, minimize,
                      accuracy, exactness);
        consistent(key.size() + 1);
      }
    } else {
      assert(!std::isinf(_cost));
      assert(_var != std::numeric_limits<uint32_t>::max());
      assert(false);
    }
    consistent(key.size() + 1);
    return;
  }

  if (prefix > key.size()) {

    if (size_t var = tree["var"].get<uint32_t>() + 1 + key.size();
        var != _var) {
      //_parent->replace(this, key, tree, action, parent, prefix);
      assert(false);
    } else {
      size_t ivar = tree["var"].get<uint32_t>();
      double bound = tree["bound"].get<double>();
      if (exactness.size() > ivar && exactness[ivar] != 0) {
        bound = std::round(bound / exactness[ivar]) * exactness[ivar];
      }
      if (_limit == bound) {
        _low->insert(key, tree["low"], action, parent, prefix + 1, minimize,
                     accuracy, exactness);
        _high->insert(key, tree["high"], action, parent, prefix + 1, minimize,
                      accuracy, exactness);
      } else {
        if (_limit > bound) {
          _low->insert(key, tree, action, parent, prefix, minimize, accuracy,
                       exactness);
          _high->insert(key, tree["high"], action, parent, prefix + 1, minimize,
                        accuracy, exactness);
        } else {
          _high->insert(key, tree, action, parent, prefix, minimize, accuracy,
                        exactness);
          _low->insert(key, tree["low"], action, parent, prefix + 1, minimize,
                       accuracy, exactness);
        }
      }
    }
    consistent(key.size() + 1);
  } else if (prefix == key.size()) {
    if (_var != prefix) {
      assert(false);
    } else {
      consistent(key.size() + 1);
      if (_limit == action) {
        _low->insert(key, tree, action, parent, prefix + 1, minimize, accuracy,
                     exactness);
        consistent(key.size() + 1);
      } else {
        assert(_low->_var > prefix);
        if (_limit > action) {
          auto tmp = std::make_shared<node_t>();
          tmp->_low = std::make_shared<node_t>();
          tmp->_cost = tmp->_low->_cost =
              (minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
          tmp->_high = shared_from_this();
          tmp->_limit = action;
          tmp->_var = prefix;
          if (_parent) {
            tmp->_parent = _parent;
            if (_parent->_high.get() == this)
              _parent->_high = tmp;
            else
              _parent->_low = tmp;
          } else {
            parent._root = tmp;
          }
          tmp->_low->_parent = tmp->_high->_parent = tmp.get();
          tmp->consistent(key.size() + 1);
          tmp->insert(key, tree, action, parent, prefix, minimize, accuracy,
                      exactness);
        } else if (_high->_var == prefix && _high->_limit < action) {
          _high->insert(key, tree, action, parent, prefix, minimize, accuracy,
                        exactness);
          consistent(key.size() + 1);
        } else if (_high->_var > prefix) {
          assert(std::isinf(_high->_cost));
          _high->_low = std::make_shared<node_t>();
          _high->_high = std::make_shared<node_t>();
          _high->_low->_cost = _high->_high->_cost =
              (minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
          _high->_limit = action;
          _high->_var = prefix;
          _high->_low->_parent = _high.get();
          _high->_high->_parent = _high.get();
          _high->_parent = this;
          _high->consistent(key.size() + 1);
          _high->insert(key, tree, action, parent, prefix, minimize, accuracy,
                        exactness);
          consistent(key.size() + 1);
        } else if (_high->_var == prefix && _high->_limit > action) {
          auto tmp = std::make_shared<node_t>();
          tmp->_low = std::make_shared<node_t>();
          tmp->_low->_cost = tmp->_cost =
              (minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
          tmp->_parent = this;
          tmp->_high = _high;
          _high = tmp;
          tmp->_low->_parent = tmp.get();
          tmp->_high->_parent = tmp.get();
          tmp->_var = prefix;
          tmp->_limit = action;
          tmp->consistent(key.size() + 1);
          tmp->insert(key, tree, action, parent, prefix, minimize, accuracy,
                      exactness);
          consistent(key.size() + 1);
        } else if (_high->_var == prefix && _high->_limit == action) {
          _high->insert(key, tree, action, parent, prefix, minimize, accuracy,
                        exactness);
        } else {
          assert(false);
        }
      }
    }
    consistent(key.size() + 1);
  } else if (prefix < key.size()) {
    if (_var != prefix) {
      //_parent->replace(this, key, tree, action, parent, prefix);
      assert(false);
    } else {
      if (_limit == key[prefix]) {
        _low->insert(key, tree, action, parent,
                     prefix + (_low->_var != _var ? 1 : 0), minimize, accuracy,
                     exactness);
      } else {
        const auto b = _limit < key[prefix];
        auto branch = (*this)[b];
        if (branch->_var != prefix) {
          // we need to inject a node here
          auto next = std::make_shared<node_t>();
          (*next)[false] = std::make_shared<node_t>();
          (*next)[false]->_cost = next->_cost =
              (minimize ? 1 : -1) * std::numeric_limits<double>::infinity();
          next->_parent = this;
          next->_var = prefix;
          next->_limit = key[prefix];
          // example
          // original [0,10] | [11,20] made for element 10,
          // we inject element 7, so we now get
          // ([0-7] | [8-10] ) | [11,20
          // where the original low is moved to the high
          // of the newly created node
          (*next)[true] = branch;
          (*this)[b] = next;
          (*next)[true]->_parent = next.get();
          (*next)[false]->_parent = next.get();
          next->insert(key, tree, action, parent, prefix, minimize, accuracy,
                       exactness);
        } else
          branch->insert(key, tree, action, parent, prefix, minimize, accuracy,
                         exactness);
      }
    }
    consistent(key.size() + 1);
  }
  consistent(key.size() + 1);
}

void SimpleTree::node_t::consistent(size_t prefix) const {
  return;
#ifndef NDEBUG
  assert(_parent != this);
  assert(_parent == nullptr || _parent->_low.get() == this ||
         _parent->_high.get() == this);
  assert(_low == nullptr || _low->_parent == this);
  assert(_high == nullptr || _high->_parent == this);
  if (_var == prefix - 1 && _parent && _parent->_var == _var) {
    assert(_parent->_limit < _limit);
  }
  if (_var < prefix)
    assert(_parent == nullptr || _parent->_var <= _var);
  if (_low)
    _low->consistent(prefix);
  if (_high)
    _high->consistent(prefix);
#endif
}

double SimpleTree::value(const double *disc, const double *cont,
                         uint32_t action) const {
  if (_root)
    return _root->value(disc, cont, action, _statevars.size());
  return std::numeric_limits<double>::infinity();
}

double SimpleTree::node_t::value(const double *disc, const double *cont,
                                 uint32_t action, size_t ndisc) const {
  if (is_leaf())
    return _cost;
  if (_var < ndisc) {
    if (disc[_var] <= _limit)
      return _low ? _low->value(disc, cont, action, ndisc) : _cost;
    return _high ? _high->value(disc, cont, action, ndisc) : _cost;
  }
  if (_var == ndisc) {
    if (action <= _limit)
      return _low ? _low->value(disc, cont, action, ndisc) : _cost;
    return _high ? _high->value(disc, cont, action, ndisc) : _cost;
  }
  if (cont[_var - (ndisc + 1)] <= _limit)
    return _low ? _low->value(disc, cont, action, ndisc) : _cost;
  return _high ? _high->value(disc, cont, action, ndisc) : _cost;
}

SimpleTree::node_ptr SimpleTree::node_t::simplify(bool make_dd,
                                                  nodemap_t &nodemap,
                                                  SimpleTree &parent) {
  if (_low)
    _low = _low->simplify(make_dd, nodemap, parent);
  if (_high)
    _high = _high->simplify(make_dd, nodemap, parent);
  if (_low)
    _low->_parent = this;
  if (_high)
    _high->_parent = this;
  if (_var < parent._statevars.size()) {
    if (_low && _low->is_leaf() && std::isinf(_low->_cost))
      return _high;
    if (_high && _high->is_leaf() && std::isinf(_high->_cost))
      return _low;
  }

  assert(_low.get() != this);
  assert(_high.get() != this);
  if ((_low == nullptr || _low->is_leaf()) &&
      (_high == nullptr || _high->is_leaf())) {
    auto lc = _low ? _low->_cost : _cost;
    auto hc = _high ? _high->_cost : _cost;
    if (lc == hc) {
      _low = nullptr;
      _high = nullptr;
      _cost = lc;
      _limit = std::numeric_limits<double>::infinity();
      _var = std::numeric_limits<decltype(_var)>::max();
      return shared_from_this();
    }
  }

  if (std::isinf(_limit) && (_low != nullptr || _high != nullptr)) {
    assert(false);
    return _low == nullptr ? _high : _low;
  }

  // if comparison on same variable in on child and other child is leaf, then
  // if the non-leaf child has a leaf-child in, we can merge.
  {
    if (_low && _low->is_leaf() && _high &&
        (!_high->is_leaf() || !std::isinf(_high->_cost)) &&
        _high->_var == _var && _high->_low && _high->_low->is_leaf() &&
        _high->_low->_cost == _low->_cost) {
      return _high;
    }

    if (_high && _high->is_leaf() && _low &&
        (!_low->is_leaf() || !std::isinf(_low->_cost)) && _low->_var == _var &&
        _low->_high && _low->_high->is_leaf() &&
        _low->_high->_cost == _high->_cost) {
      return _low;
    }
  }

  /*if(make_dd)
  {
      if(_high == nullptr)
          return _low;
      if(_low == nullptr)
          return _high;
      if(_low && _parent && !_low->is_leaf() && _high->is_leaf())
      {
          auto cost = _high ? _high->_cost : _cost;
          auto lcost = _low->_high ? _low->_high->_cost : _low->_cost;
          if(cost == lcost && _low->_var == _var)
              return _low;
      }
      if(_high && _parent && !_high->is_leaf() && _low->is_leaf())
      {
          auto cost = _low ? _low->_cost : _cost;
          auto hcost = _high->_low ? _high->_low->_cost : _high->_cost;
          if(cost == hcost && _high->_var == _var)
              return _high;
      }
  }*/
  if (_high != nullptr && _high == _low) {
    return _high;
  }
  if (make_dd) {
    auto &ptr = nodemap[*this];
    if (ptr == nullptr)
      ptr = shared_from_this();
    return ptr;
  }
  return shared_from_this();
}

/*void SimpleTree::node_t::rec_insert(double value,
std::vector<std::pair<bool,bool>>& handled, std::vector<std::pair<double,
double> >& bounds, bool minimize) { auto nid = _var; assert(nid <
bounds.size()); bool reset_low = false; bool reset_high = false; bool best =
false; if( std::isnan(_cost)            || (minimize && value < _cost)  ||
        (!minimize && value > _cost))
    {
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
            _low->_limit = (!handled[nid].first ? bounds[nid].first :
bounds[nid].second);
        }
        _low->rec_insert(value, handled, bounds, minimize);
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
            _high->_limit = (!handled[nid].first ? bounds[nid].first :
bounds[nid].second);
        }
        _high->rec_insert(value, handled, bounds, minimize);
        if(reset_low) handled[_var].first = false;
        if(reset_high) handled[_var].second = false;
    }
    else
    {
        if(best && handled[_var].first && handled[_var].second && _var ==
handled.size()-1) { _low = nullptr; _high = nullptr; if(reset_low)
handled[_var].first = false; if(reset_high) handled[_var].second = false;
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
            _low->_limit = (!handled[nid].first ? bounds[nid].first :
bounds[nid].second);
        }
        if(_high == nullptr)
        {
            _high = std::make_shared<node_t>();
            _high->_parent = this;
            _high->_var = nid;
            _high->_limit = (!handled[nid].first ? bounds[nid].first :
bounds[nid].second);
        }

        _low->rec_insert(value, handled, bounds, minimize);
        _high->rec_insert(value, handled, bounds, minimize);
    }
    if(reset_low) handled[_var].first = false;
    if(reset_high) handled[_var].second = false;
}*/

std::ostream &SimpleTree::print(std::ostream &stream) const {
  _root->print(stream, 0);
  return stream;
}

std::ostream &SimpleTree::node_t::print(std::ostream &out, size_t tabs) const {
  for (size_t i = 0; i < tabs; ++i)
    out << "\t";
  out << "{\"var\":" << _var << ",\"bound\":" << _limit;
  if (!std::isnan(_cost)) {
    out << ",\"value\":" << _cost;
  }
  if (_low) {
    out << ",\n";
    for (size_t i = 0; i < tabs; ++i)
      out << "\t";
    out << "\"low\":\n";
    _low->print(out, tabs + 1);
  }
  if (_high) {
    out << ",\n";
    for (size_t i = 0; i < tabs; ++i)
      out << "\t";
    out << "\"high\":\n";
    _high->print(out, tabs + 1);
  }
  out << "\n";
  for (size_t i = 0; i < tabs; ++i)
    out << "\t";
  out << "}";
  return out;
}

std::ostream &SimpleTree::print_c(std::ostream &os,
                                  const std::string &name) const {
  if (_root == nullptr) {
    os << "double " << name
       << "(unsigned int action, const double* disc, const double* vars)\n{\n";
    os << "\treturn 0 ; \n}\n";
    return os;
  }
  auto _idnode = std::unordered_map<size_t, const node_t *>{};
  auto _nodeid = std::unordered_map<const node_t *, size_t>{};
  _idnode[0] = _root.get();
  _nodeid[_root.get()] = 0;

  size_t lid = 0;
  while (lid != _nodeid.size()) {
    auto n = _idnode[lid];
    if (!n->is_leaf()) {
      auto ln = n->_low.get();
      auto hn = n->_high.get();
      if (_nodeid.count(ln) == 0) {
        auto id = _idnode.size();
        _idnode[id] = ln;
        _nodeid[ln] = id;
      }
      if (_nodeid.count(hn) == 0) {
        auto id = _idnode.size();
        _idnode[id] = hn;
        _nodeid[hn] = id;
      }
    }
    ++lid;
  }

  os << "const int " << name << "_nodes[] = {";
  for (size_t i = 0; i < lid; ++i) {
    if (i != 0)
      os << ",";
    auto n = _idnode[i];
    if (n->is_leaf()) {
      os << "-1,-1,-1";
    } else {
      os << _nodeid[n->_low.get()] << ",";
      os << _nodeid[n->_high.get()] << ",";
      os << n->_var;
    }
  }
  os << "};\n";
  const auto [mmin, mmax] = _root->compute_min_max();
  auto v = is_minimization() ? mmax + 1 : mmin - 1;
  os << "const double " << name << "_values[] = {";
  for (size_t i = 0; i < lid; ++i) {
    if (i != 0)
      os << ",";
    auto n = _idnode[i];
    if (n->is_leaf()) {
      if (!std::isinf(n->_cost) && !std::isnan(n->_cost))
        os << n->_cost;
      else
        os << v;
    } else {
      os << n->_limit;
    }
  }
  os << "};\n";
  os << "double " << name
     << "(unsigned int action, const double* disc, const double* vars)\n{\n";
  // stream << "\t// Depth = " << _root->depth() << std::endl;
  os << "\t// Actions = " << _actions.size() << std::endl;
  os << "\t// Disc = " << _statevars.size() << std::endl;
  os << "\t// Cont = " << _pointvars.size() << std::endl;
  os << "\t// Nodes = " << lid << std::endl;
  /*std::unordered_set<const node_t*> printed;
  std::vector<const node_t*> toprint;
  if(_root)
      _root->print_c_nested(stream, _statevars.size(), 1, toprint, _root);
  auto mm = _root->compute_min_max();
  auto v = is_minimization() ? mm.second + 1 : mm.first -1;
  stream << "\treturn " << v << ";\n";
  for(auto n : toprint)
  {
      n->print_c(stream, 0, printed, 1);
      stream << "\treturn " << v << ";\n";
  }
  stream << "\treturn " << v << ";\n";
   * */
  os << "\tint ins = 0;\n\twhile(true) {\n";
  os << "\t\tint l = " << name << "_nodes[ins*3]; int h = " << name
     << "_nodes[1+(ins*3)]; int v = " << name << "_nodes[2+(ins*3)];\n";
  os << "\t\tif(v == -1) return " << name << "_values[ins];\n";
  os << "\t\tdouble val = 0;\n";
  os << "\t\tif(v == " << _statevars.size() << ") val = action;\n";
  os << "\t\telse if(v > " << _statevars.size() << ") val = vars[v-"
     << (_statevars.size() + 1) << "];\n";
  os << "\t\telse val = disc[v];\n";
  os << "\t\tif(val <= " << name << "_values[ins])\n";
  os << "\t\t\tins = l;\n";
  os << "\t\telse\n";
  os << "\t\t\tins = h;\n";
  os << "\t}\n";
  os << "\treturn " << v << ";\n";
  os << "}\n";
  return os;
}

size_t SimpleTree::node_t::depth() const {
  if (_low && _high == nullptr)
    return 1 + _low->depth();
  if (_high && _low == nullptr)
    return 1 + _high->depth();
  if (_low == nullptr && _high == nullptr)
    return 0;
  return 1 + std::max(_low->depth(), _high->depth());
}

std::ostream &
SimpleTree::node_t::print_c(std::ostream &os, size_t disc,
                            std::unordered_set<const node_t *> &printed,
                            size_t tabs) const {

  if (printed.count(this) > 0)
    return os;
  printed.insert(this);
  // for(size_t i = 0; i < tabs; ++i) stream << "\t";
  os << "l" << this << ":\n";
  std::vector<const node_t *> toprint;
  if (is_leaf()) {
    //    for(size_t i = 0; i < tabs+1; ++i) stream << "\t";
    if (!std::isinf(_cost))
      os << "return " << _cost << ";" << std::endl;
    else
      os << "{}";
    return os;
  }
  if (std::isinf(_limit)) {
    if (_low)
      return _low->print_c(os, disc, printed, tabs);
    if (_high)
      return _high->print_c(os, disc, printed, tabs);
    return os;
  }
  // for(size_t i = 0; i < tabs+1; ++i) stream << "\t";
  if (_var == disc) {
    os << "if(action <= " << _limit << ") ";
  } else {
    auto type = disc > _var ? "disc" : "vars";
    auto vid = disc > _var ? _var : _var - (disc + 1);
    os << "if(" << type << "[" << vid << "] <= " << _limit << ") ";
  }
  if (_low)
    _low->print_c_nested(os, disc, tabs + 1, toprint, _low);
  else {
    if (!std::isinf(_cost))
      os << "return " << _cost << ";";
    else
      os << "{}";
  }
  os << "\n";
  //    for(size_t i = 0; i < tabs+1; ++i) stream << "\t";
  os << "else ";
  if (_high)
    _high->print_c_nested(os, disc, tabs + 1, toprint, _high);
  else {
    if (!std::isinf(_cost))
      os << "return " << _cost << ";";
    else
      os << "{}";
  }
  os << "\n";
  for (auto n : toprint)
    n->print_c(os, disc, printed, tabs);
  return os;
}

std::ostream &
SimpleTree::node_t::print_c_nested(std::ostream &os, size_t disc, size_t tabs,
                                   std::vector<const node_t *> &toprint,
                                   const node_ptr &node) const {
  assert(node.get() == this);
  if (is_leaf() && std::isinf(_limit)) {
    if (!std::isinf(_cost))
      os << "return " << _cost << ";";
    else
      os << "{}";
    return os;
  }
  if (node.use_count() == 1) {
    os << "{\n";
    //        for(size_t i = 0; i < tabs; ++i) stream << "\t";
    if (_var == disc) {
      os << "if(action <= " << _limit << ") ";
    } else {
      auto type = disc > _var ? "disc" : "vars";
      auto vid = disc > _var ? _var : _var - (disc + 1);
      os << "if(" << type << "[" << vid << "] <= " << _limit << ") ";
    }
    if (_low)
      _low->print_c_nested(os, disc, tabs + 1, toprint, _low);
    else {
      if (!std::isinf(_cost))
        os << "return " << _cost << ";";
      else
        os << "{}";
    }
    os << "\n";
    // for(size_t i = 0; i < tabs; ++i) stream << "\t";
    os << "else ";
    if (_high)
      _high->print_c_nested(os, disc, tabs + 1, toprint, _high);
    else {
      if (!std::isinf(_cost))
        os << "return " << _cost << ";";
      else
        os << "{}";
    }
    os << "\n";
    // for(size_t i = 0; i < tabs-1; ++i) stream << "\t";
    os << "}";
    return os;
  }
  os << "goto l" << this << ";";
  toprint.push_back(this);
  return os;
}
