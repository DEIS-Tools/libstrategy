#include <libstrategy/SimpleTree.h>

#include <fstream>
#include <iostream>

int main(int argc, const char* argv[])
{
    if (argc < 2) {
        std::cout << "Expecting a path to a strategy as an argument\n";
        return 1;
    }
    for (auto i = 1; i < argc; ++i) {
        auto strategy_path = argv[i];
        auto is = std::ifstream{strategy_path};
        auto tree = SimpleTree::parse(is, true, false);
        tree.print(std::cout);
    }
}
