#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <charconv>    // from_chars
#include <type_traits> // true_type for detecting from_chars
#include <utility>     // declval for detecting from_chars

/// C++17 compile-time test for presence of std::from_chars(const char*, const
/// char*, T&) Replace it with C++20 concepts later (or perhaps AppleClang will
/// implement proper from_chars by then). History: C++17 introduced
/// std::from_chars, but STL vendors were late, then provided only integral
/// versions...
template <typename, typename = void>
struct has_from_chars : std::false_type {
}; // primary template declaration (used when specializations fail)

template <typename T>
struct has_from_chars< // template partial specialization
    T,
    std::void_t< // tests if the following expression computes into a type:
        decltype(std::from_chars(std::declval<const char *&>(),
                                 std::declval<const char *&>(),
                                 std::declval<T &>()))>> : std::true_type {};

template <typename T>
constexpr auto has_from_chars_v = has_from_chars<T>::value;

#endif // UTILITIES_HPP
