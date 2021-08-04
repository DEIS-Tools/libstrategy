/*
 * Copyright (C) 2021 Peter G. Jensen <root@petergjoel.dk>
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

#define BOOST_TEST_MODULE UnorderedLoad

#include <boost/test/unit_test.hpp>
#include <fstream>

#include "SimpleTree.h"

BOOST_AUTO_TEST_CASE(DirectoryTest)
{
    BOOST_REQUIRE(getenv("STRATEGY_DIR"));
}

BOOST_AUTO_TEST_CASE(Inconsistent1)
{
    std::string strategy = getenv("STRATEGY_DIR");
    strategy += "/inconsistent1.strategy";
    std::ifstream in(strategy);
    auto tree = SimpleTree::parse(in, false, false);
    double vars[] = {10};
    auto act18 = tree.value(vars,nullptr, 0);
    auto act19 = tree.value(vars,nullptr, 1);
    BOOST_REQUIRE_LT(act18, act19);
}

BOOST_AUTO_TEST_CASE(Inconsistent1Simplify)
{
    std::string strategy = getenv("STRATEGY_DIR");
    strategy += "/inconsistent1.strategy";
    std::ifstream in(strategy);
    auto tree = SimpleTree::parse(in, true, false);
    double vars[] = {10};
    BOOST_REQUIRE_LT(tree.value(vars,nullptr, 0), tree.value(vars,nullptr, 1));
}

BOOST_AUTO_TEST_CASE(Inconsistent1SimplifySubsumption)
{
    std::string strategy = getenv("STRATEGY_DIR");
    strategy += "/inconsistent1.strategy";
    std::ifstream in(strategy);
    auto tree = SimpleTree::parse(in, true, true);
    double vars[] = {10};
    BOOST_REQUIRE_LT(tree.value(vars,nullptr, 0), tree.value(vars,nullptr, 1));
}
