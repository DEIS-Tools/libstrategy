#define BOOST_TEST_MODULE UnorderedLoad

#include "SimpleTree.h"
#include "errors.h"
#include "utilities.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(SimpleTreeParseKey) {
  const auto res_blank = parse_key("()");
  BOOST_CHECK(res_blank.empty());

  const auto res_one = parse_key("(1)");
  BOOST_REQUIRE_EQUAL(res_one.size(), 1);
  BOOST_CHECK_EQUAL(res_one[0], 1);

  const auto res_two = parse_key("(2,1)");
  BOOST_REQUIRE_EQUAL(res_two.size(), 2);
  BOOST_CHECK_EQUAL(res_two[0], 2);
  BOOST_CHECK_EQUAL(res_two[1], 1);

  const auto res_three = parse_key("(3,2,1)");
  BOOST_REQUIRE_EQUAL(res_three.size(), 3);
  BOOST_CHECK_EQUAL(res_three[0], 3);
  BOOST_CHECK_EQUAL(res_three[1], 2);
  BOOST_CHECK_EQUAL(res_three[2], 1);

  const auto res_float = parse_key("(3.141)");
  BOOST_REQUIRE_EQUAL(res_float.size(), 1);
  BOOST_CHECK_EQUAL(res_float[0], 3.141);

  const auto res_floats = parse_key("(4.3,2.1)");
  BOOST_REQUIRE_EQUAL(res_floats.size(), 2);
  BOOST_CHECK_EQUAL(res_floats[0], 4.3);
  BOOST_CHECK_EQUAL(res_floats[1], 2.1);

  // a few negative tests:
  BOOST_CHECK_THROW(parse_key(""), base_error);
  BOOST_CHECK_THROW(parse_key("1"), base_error);
  BOOST_CHECK_THROW(parse_key(")"), base_error);
  BOOST_CHECK_THROW(parse_key("("), base_error);
  BOOST_CHECK_THROW(parse_key("(1"), base_error);
  BOOST_CHECK_THROW(parse_key("(2,"), base_error);
  BOOST_CHECK_THROW(parse_key("(2,1"), base_error);
}
