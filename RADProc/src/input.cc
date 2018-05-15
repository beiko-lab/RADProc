// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

//
// input.cc -- routines to read various formats of data into the XXX data structure.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//

#include "input.h"

Seq::Seq() { 
    this->id       = NULL;
    this->seq      = NULL;
    this->qual     = NULL;
    this->loc_str  = NULL;
}


Input::Input() {
    memset(this->line, '\0', max_len);
}

Input::Input(const char *path) {
    memset(this->line, '\0', max_len);

    this->path = string(path);
    //
    // Open the file for reading
    //
    this->fh.open(path, ifstream::in);

    if (this->fh.fail()) 
        cerr << "Error opening input file '" << path << "'\n";
}

Input::~Input() {
    // Close the file
    this->fh.close();
}


