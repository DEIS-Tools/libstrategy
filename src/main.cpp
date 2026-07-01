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
 * Author: Peter G. Jensen
 *
 * Created on December 13, 2018, 2:11 PM
 */

#include "ZonotopStrategy.h"

#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <fstream>
#include <cctype>
#include <cstdlib>
#include <charconv>

namespace {

void print_help(const char* prog)
{
    std::cout << "Usage: " << prog
              << " [options]\n\n"
                 "Options:\n"
                 "  -h, --help                     produce help message\n"
                 "  -i, --input <file>             Input of synthesized controller.\n"
                 "  -l, --learned <file>           Input of learned controller.\n"
                 "  -c, --choice <name>            Output function name for choice.\n"
                 "  -p, --pattern <name>           Output function name for pattern.\n"
                 "  -a, --accuracy <value>\n"
                 "  -e, --exactness <value>...\n";
}

// A token counts as a value (not the start of a new option) unless it
// looks like "-x" or "--xxx" (as opposed to a negative number like "-1.5").
bool is_option(std::string_view s)
{
    return s.size() > 1 && s[0] == '-' && !(std::isdigit(static_cast<unsigned char>(s[1])) || s[1] == '.');
}

double to_double(std::string_view s)
{
    double value;
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
    if (ec != std::errc() || ptr != s.data() + s.size()) {
        std::cerr << "Invalid number: " << s << std::endl;
        std::exit(-1);
    }
    return value;
}

}  // namespace

int main(int argc, char** argv)
{
    std::string infile;
    std::string learned;
    std::string learned_name = "choice";
    std::string pattern_name = "pattern";
    double accuracy = 0;
    std::vector<double> exactness;

    std::vector<std::string_view> args(argv + 1, argv + argc);
    for (size_t i = 0; i < args.size(); ++i) {
        std::string_view arg = args[i];
        auto next_value = [&]() -> std::string_view {
            if (i + 1 >= args.size() || is_option(args[i + 1])) {
                std::cerr << "Missing value for option " << arg << std::endl;
                std::exit(-1);
            }
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") {
            print_help(argv[0]);
            return 1;
        } else if (arg == "-i" || arg == "--input") {
            infile = next_value();
        } else if (arg == "-l" || arg == "--learned") {
            learned = next_value();
        } else if (arg == "-c" || arg == "--choice") {
            learned_name = next_value();
        } else if (arg == "-p" || arg == "--pattern") {
            pattern_name = next_value();
        } else if (arg == "-a" || arg == "--accuracy") {
            accuracy = to_double(next_value());
        } else if (arg == "-e" || arg == "--exactness") {
            while (i + 1 < args.size() && !is_option(args[i + 1])) {
                exactness.push_back(to_double(args[++i]));
            }
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_help(argv[0]);
            return -1;
        }
    }
    std::cout << "#include \"sub_ctrl.h\"\n";
    if (!infile.empty()) {
        std::ifstream instream(infile);
        if (instream.fail()) {
            std::cerr << "Could not open infile for reading : " << infile << std::endl;
            return -1;
        }
        auto strategy = ZonotopStrategy::parse(instream);
        std::cout << "// From \"" << infile << "\"" << std::endl;
        strategy.print_c(std::cout, pattern_name);
    }
    if (learned.empty()) {
        return 0;
    }
    std::ifstream lstream(learned);
    if (lstream.fail()) {
        std::cerr << "Could not open learned for reading : " << learned << std::endl;
        return -1;
    }
    std::cout << "// From \"" << learned << "\"" << std::endl;
    std::cerr << "EXACTNESS ";
    for (auto e : exactness)
        std::cerr << e << " , ";
    std::cerr << std::endl;
    auto learned_strategy = SimpleTree::parse(lstream, true, true, accuracy, exactness);
    // strategy.filter(learned_strategy);
    learned_strategy.print_c(std::cout, learned_name);
    return 0;
}
