/*
 *  Copyright Peter G. Jensen, all rights reserved.
 */

/* 
 * File:   libz2s.h
 * Author: Peter G. Jensen <root@petergjoel.dk>
 *
 * Created on May 8, 2019, 10:34 AM
 */

#ifndef LIBZ2S_H
#define LIBZ2S_H
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif

void* parse_strategy(const char*);
void destroy_strategy(void*);

int num_patterns(void*);
int max_pattern_length(void*);
void active(void*, const double* sample, bool* write);
int get_pattern(void*, int el, int* write);
double get_min(void*, int dimen);
double get_max(void*, int dimen);

#ifdef __cplusplus
}
#endif

#endif /* LIBZ2S_H */

