// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __DNASeq_H__
#define __DNASeq_H__

#include <string.h>
#include <limits.h>

// #ifdef __GNUC__
// #include <ext/hash_map>
// using __gnu_cxx::hash_map;
// using __gnu_cxx::hash;
// #else
// #include <hash_map>
// #endif

//
// We expect (and C++ defines) an unsigned char as 8 bits, so we 
// should be able to store 4 nucleotide bases per byte of memory.
//
const unsigned short int bases_per_byte = CHAR_BIT / 2;

//
// DNA Sequence Storage Class
//
// Two-bit compression, four bases per byte of storage:
//    A == 00
//    C == 01
//    G == 10
//    T == 11
//
class DNASeq {
public:
    //
    // The number of DNA bases we are storing
    //
    unsigned short int size;
    //
    // Array of bytes to store DNA sequence, one character per two bits, four per byte.
    //
    unsigned char *s;

    DNASeq(int);
    DNASeq(int, const char *);
    DNASeq(int, unsigned char *);
    ~DNASeq();

    char  len() { return this->size; }
    char  operator[](int);
    char *seq(char *);
    char *seq();
    char *subseq(char *, int, int);
};

#include <iostream>
#include <fstream>
#include <sstream>
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;

struct hash_dnaseq {
    size_t operator()(DNASeq *__s) const
    {
	size_t __result = static_cast<size_t>(14695981039346656037ULL);
	unsigned short int   __bytes  = (__s->size / bases_per_byte) + (__s->size % bases_per_byte > 0 ? 1 : 0);
	for (unsigned short int i = 0; i < __bytes; i++) {
	    __result ^= static_cast<size_t>(__s->s[i]);
	    __result *= static_cast<size_t>(1099511628211ULL);
	}

	return __result;
    }
};

struct dnaseq_eqstr {
    bool operator()(DNASeq *s1, DNASeq *s2) const {
	unsigned int bytes = (s1->size / bases_per_byte) + (s1->size % bases_per_byte > 0 ? 1 : 0);
	for (unsigned int i = 0; i < bytes; i++)
	    if (s1->s[i] != s2->s[i]) return false;
	return true;
    }
};

#endif // __DNASeq_H__
