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

#ifndef __INPUT_H__
#define __INPUT_H__

#include <errno.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using std::ifstream;
using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include "constants.h"



class Seq {
 public:
    char *id;
    char *seq;
    char *qual;

    // Location information for a mapped sequence
    char  *loc_str;
    //PhyLoc loc;

    Seq( void );
    
};

//
// The base class for all of our Input classes, such as Tsv, Fastq, Fasta, etc.
//
class Input {
 public:
    string   path;
    ifstream fh;
    char     line[max_len];

    Input();
    Input(const char *path);
    virtual ~Input();
    virtual Seq *next_seq() = 0;
    virtual int  next_seq(Seq &) = 0;
};


#endif // __INPUT_H__
