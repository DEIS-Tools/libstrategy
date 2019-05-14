/* 
 * Author: Peter G. Jensen
 *
 * Created on December 13, 2018, 2:11 PM
 */

#ifndef BOUND_H
#define BOUND_H

#include <limits>


struct bound_t {
    double _min = std::numeric_limits<double>::infinity();
    double _max = -std::numeric_limits<double>::infinity();
};

#endif /* BOUND_H */

