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

void* parse_learned(const char*);
void destroy_learned(void*);
double weight(void*, const double* disc, const double* cont, int a);

#ifdef __cplusplus
}
#endif

#endif /* LIBZ2S_H */

