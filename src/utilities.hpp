#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "errors.h"

/**
 * Avoid including this header into headers (include into cpp instead)
 * as it includes streams (instead of iosfwd) and slows down compilation.
 */

#include <charconv>     // from_chars
#include <sstream>      // fallback if from_chars is not available
#include <type_traits>  // true_type for detecting from_chars
#include <utility>      // declval for detecting from_chars
#include <vector>

/** C++17 compile-time test for presence of std::from_chars(const char*, const char*, T&)
 * Replace it with C++20 concepts later (or perhaps AppleClang will implement proper from_chars by then).
 * History: C++17 introduced std::from_chars, but STL vendors were late, then provided only integral versions...
 */
template <typename, typename = void>
struct has_from_chars : std::false_type
{};  // primary template declaration (used when specializations fail)

template <typename T>
struct has_from_chars<  // template partial specialization
    T,
    std::void_t<  // tests if the following expression computes into a type:
        decltype(std::from_chars(std::declval<const char*&>(), std::declval<const char*&>(), std::declval<T&>()))>>
    : std::true_type
{};

template <typename T>
constexpr auto has_from_chars_v = has_from_chars<T>::value;

/// Parses numbers from a key in a form of "(number,number,number)"
template <typename T = double>  // has to be a template, otherwise AppleClang ignores constexpr
std::vector<T> parse_key(const std::string& key)
{
    static_assert(std::is_arithmetic_v<T>, "only numeric keys are supported");
    auto res = std::vector<T>{};
    T number;                             // the number to parse into
    if constexpr (has_from_chars_v<T>) {  // fast floating point parsing
        auto it = key.c_str();
        const auto end = it + key.size();
        if (it == end || *it != '(')
            throw base_error("incorrectly formatted key ('(' expected): " + key);
        if (*++it == ')')
            return res;
        while (it != end) {
            if (auto [p, ec] = std::from_chars(it, end, number); ec == std::errc()) {
                res.push_back(number);
                it = p;
                if (it != end && *it == ',')
                    ++it;
                else
                    break;
            } else
                throw base_error("failed to parse number in key: " + key);
        }
        if (it == end || *it != ')')
            throw base_error("incorrectly formatted key (')' expected): " + key);
    } else {  // fallback to slow stream parsing (AppleClang does not support from_chars)
        auto is = std::istringstream{key};
        char c;
        if (!is.get(c) || c != '(')
            throw base_error("incorrectly formatted key ('(' expected): " + key);
        if (is && is.peek() == ')')
            return res;
        while (is >> number) {
            res.push_back(number);
            if (!is.get(c) || c != ',')
                break;
        }
        if (c != ')') {
            throw base_error("incorrectly formatted key (')' expected): " + key);
        }
    }
    return res;
}

#endif  // UTILITIES_HPP
