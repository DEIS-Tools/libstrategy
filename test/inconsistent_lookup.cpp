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

#include "SimpleTree.h"

#include <doctest/doctest.h>

#include <fstream>
#include <filesystem>

TEST_SUITE_BEGIN("Unordered Load");

namespace fs = std::filesystem;

TEST_CASE("Inconsistent1")
{
    const auto strategy_dir = getenv("STRATEGY_DIR");
    REQUIRE(strategy_dir != nullptr);

    const auto strategy_path = fs::path{strategy_dir} / "inconsistent1.strategy";
    auto in = std::ifstream{strategy_path};
    double vars[] = {10};
    SUBCASE("Plain")
    {
        auto tree = SimpleTree::parse(in, false, false);
        auto act18 = tree.value(vars, nullptr, 0);
        auto act19 = tree.value(vars, nullptr, 1);
        REQUIRE(act18 < act19);
    }
    SUBCASE("Simplify")
    {
        auto tree = SimpleTree::parse(in, true, false);
        REQUIRE(tree.value(vars, nullptr, 0) < tree.value(vars, nullptr, 1));
    }
    SUBCASE("Simplify Subsumption")
    {
        auto tree = SimpleTree::parse(in, true, true);
        REQUIRE(tree.value(vars, nullptr, 0) < tree.value(vars, nullptr, 1));
    }
}

TEST_SUITE_END();
