#include "SimpleTree.h"
#include "errors.h"
#include "utilities.hpp"

#include <doctest/doctest.h>

TEST_CASE("SimpleTree::parse_key")
{
    const auto res_blank = parse_key("()");
    CHECK(res_blank.empty());

    const auto res_one = parse_key("(1)");
    REQUIRE(res_one.size() == 1);
    CHECK(res_one[0] == 1);

    const auto res_two = parse_key("(2,1)");
    REQUIRE(res_two.size() == 2);
    CHECK(res_two[0] == 2);
    CHECK(res_two[1] == 1);

    const auto res_three = parse_key("(3,2,1)");
    REQUIRE(res_three.size() == 3);
    CHECK(res_three[0] == 3);
    CHECK(res_three[1] == 2);
    CHECK(res_three[2] == 1);

    const auto res_float = parse_key("(3.141)");
    REQUIRE(res_float.size() == 1);
    CHECK(res_float[0] == 3.141);

    const auto res_floats = parse_key("(4.3,2.1)");
    REQUIRE(res_floats.size() == 2);
    CHECK(res_floats[0] == 4.3);
    CHECK(res_floats[1] == 2.1);

    // a few negative tests:
    CHECK_THROWS_AS(parse_key(""), base_error);
    CHECK_THROWS_AS(parse_key("1"), base_error);
    CHECK_THROWS_AS(parse_key(")"), base_error);
    CHECK_THROWS_AS(parse_key("("), base_error);
    CHECK_THROWS_AS(parse_key("(1"), base_error);
    CHECK_THROWS_AS(parse_key("(2,"), base_error);
    CHECK_THROWS_AS(parse_key("(2,1"), base_error);
}
