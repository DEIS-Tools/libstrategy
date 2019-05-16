/* 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Author: Peter G. Jensen
 *
 * Created on December 13, 2018, 2:11 PM
 */

#include "ZonotopStrategy.h"

#include <boost/program_options.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

namespace po = boost::program_options;

int main(int argc, char **argv) {

    po::options_description opts;
    std::string infile;
    std::string learned;
    std::string learned_name = "choice";
    std::string pattern_name = "pattern";
    opts.add_options()
            ("help,h", "produce help message")
            ("input,i", po::value<std::string>(&infile),
                "Input of synthesized controller.")
            ("learned,l", po::value<std::string>(&learned),
                "Input of learned controller.")
            ("choice,c", po::value<std::string>(&learned_name),
                "Output function name for choice.")
            ("pattern,p", po::value<std::string>(&pattern_name),
                "Output function name for pattern.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << "\n";
        return 1;
    }
    size_t split = std::numeric_limits<size_t>::max();
    if(!infile.empty())
    {
        std::ifstream instream(infile);
        if(instream.fail())
        {
            std::cerr << "Could not open infile for reading : " << infile << std::endl;
            return -1;            
        }
        auto strategy = ZonotopStrategy::parse(instream);
        strategy.print_c(std::cout, pattern_name);
        split = strategy.num_patterns();
    }    
    if(learned.empty())
    {
        return 0;
    }
    std::ifstream lstream(learned);
    if(lstream.fail())
    {
        std::cerr << "Could not open learned for reading : " << learned << std::endl;
        return -1;            
    }
    auto learned_strategy = SimpleTree::parse(lstream);
    //strategy.filter(learned_strategy);
    learned_strategy.print_c(std::cout, learned_name, split);
    return 0;
}

