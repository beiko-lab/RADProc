// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2018, Praveen Nadukkalam Ravindran <pravindran@dal.ca>
//
// This file is part of RADProc.
//
// RADProc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RADProc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RADProc.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef __NETWORK_H__
#define __NETWORK_H__

#include "constants.h" 

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif

#include <getopt.h> // Process command-line options
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;

#include <iostream>
#include <fstream>
#include <sstream>
using std::ofstream;
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <unordered_map>
using std::unordered_map;
#include <queue>
using std::queue;
#include <set>
using std::set;
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif

#include "tags.h"
#include "DNASeq.h"     // Class for storing two-bit compressed DNA sequences
#include "FastaI.h"     // Reading input files in FASTA format
#include "FastqI.h"     // Reading input files in FASTQ format
#include "gzFasta.h"    // Reading gzipped input files in FASTA format
#include "gzFastq.h"    // Reading gzipped input files in FASTQ format

typedef unsigned int uint;

class HVal {
 public:
    vector<int> ids;

    int count() {
	return this->ids.size();
    }
    int add_id(int id) {
    	this->ids.push_back(id);
	return 0;
    }
};

struct hash_charptr {
    size_t operator()(const char *__s) const
    {
        size_t __result = static_cast<size_t>(14695981039346656037ULL);
        unsigned int __len = strlen(__s);
        for (unsigned int i = 0; i < __len; i++) {
            __result ^= static_cast<size_t>(__s[i]);
            __result *= static_cast<size_t>(1099511628211ULL);
        }

        return __result;
    }
};
struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
        return strcmp(s1, s2) == 0;
    }
};

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<DNASeq *, HVal, hash_dnaseq, dnaseq_eqstr> DNASeqHashMap;
typedef sparse_hash_map<const char *, vector<string>, hash_charptr, eqstr> KmerHashMap;
#else
typedef unordered_map<DNASeq *, HVal, hash_dnaseq, dnaseq_eqstr> DNASeqHashMap;
typedef unordered_map<const char *, vector<string>, hash_charptr, eqstr> KmerHashMap;
#endif

//int  load_radtags(DNASeqHashMap &, vector<DNASeq *> &, string);
int  load_radtags(map<string,int> &, string);



#endif // __USTACKS_H__
