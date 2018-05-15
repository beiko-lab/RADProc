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
// DNASeq.cc
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: DNASeq.cc 2133 2011-06-07 04:07:41Z catchen $
//

#include "DNASeq.h"

DNASeq::DNASeq(int size) {
    int bytes;

    this->size = size;

    bytes  = size / bases_per_byte;
    bytes += size % bases_per_byte > 0 ? 1 : 0;

    this->s = new unsigned char[bytes];
    memset(this->s, 0, bytes);
}

DNASeq::DNASeq(int size, unsigned char *seq) {
    unsigned int bytes;

    this->size = size;

    bytes  = size / bases_per_byte;
    bytes += size % bases_per_byte > 0 ? 1 : 0;

    this->s = new unsigned char[bytes];
    for (unsigned int i = 0; i < bytes; i++)
	this->s[i] = seq[i];
}

DNASeq::DNASeq(int size, const char *seq) {
    int bytes, rem;

    this->size = size;

    bytes  = size / bases_per_byte;
    rem    = size % bases_per_byte;
    bytes += rem > 0 ? 1 : 0;

    this->s = new unsigned char[bytes];

    int index = 0;

    for (int i = 0; i < this->size; i++) {
	//cerr << "Encoding character " << i << ", '" << seq[i] << "'\n";

	if (i > 0 && i % bases_per_byte == 0) index++;

	//cerr << "  encoding '" << seq[i] << "' into byte " << index << ".\n";

	this->s[index] <<= 2;

	switch (seq[i]) {
	case 'A':
	case 'a':
	    // A == 00
	    break;
	case 'C':
	case 'c':
	    // C == 01
	    this->s[index] |= 0x1;
	    break;
	case 'G':
	case 'g':
	    // G == 10
	    this->s[index] |= 0x2;
	    break;
	case 'T':
	case 't':
	    // T == 11
	    this->s[index] |= 0x3;
	    break;
	}

	//cerr << "    s[" << index << "," << i % bases_per_byte << "] == " << (int)this->s[index] << "\n";
    }

    if (rem > 0)
	this->s[index] <<= (bases_per_byte - rem) * 2;
}

DNASeq::~DNASeq() {
    delete [] this->s;
}

char DNASeq::operator[](int pos) {
    unsigned char c, base;
    int index, rem;

    if (pos > (this->size - 1)) return '\0';

    index = pos / bases_per_byte;
    rem   = pos % bases_per_byte;

    //cerr << "s[" << index << "," << rem << "] == " << (int)this->s[index] << "\n";

    switch (rem) {
    case 0:
	c = this->s[index] & 0xC0; // 11000000
	c >>= 6;
	break;
    case 1:
	c = this->s[index] & 0x30; // 00110000
	c >>= 4;
	break;
    case 2:
	c = this->s[index] & 0xC;  // 00001100
	c >>= 2;
	break;
    case 3:
	c = this->s[index] & 0x3;  // 00000011
	break;
    }

    switch (c) {
    case 0:
	base = 'A';
	break;
    case 1:
	base = 'C';
	break;
    case 2:
	base = 'G';
	break;
    case 3:
	base = 'T';
	break;
    }
    //cerr << "  Decoding character " << pos << ", '" << base << "'\n";

    return base;
}

char *DNASeq::subseq(char *seq, int start, int end) {
    unsigned char c;
    int i, index, rem;

    i      = start;
    index  = i / bases_per_byte;
    rem    = i % bases_per_byte;

    for (; i <= end; i++) {
	rem = i % bases_per_byte;

	if (i > 0 && rem == 0) index++;
	//cerr << "s[" << index << "," << rem << "] == " << (int)this->s[index] << "\n";

	switch (rem) {
	case 0:
	    c = this->s[index] & 0xC0; // 11000000
	    c >>= 6;
	    break;
	case 1:
	    c = this->s[index] & 0x30; // 00110000
	    c >>= 4;
	    break;
	case 2:
	    c = this->s[index] & 0xC;  // 00001100
	    c >>= 2;
	    break;
	case 3:
	    c = this->s[index] & 0x3;  // 00000011
	    break;
	}

	switch (c) {
	case 0:
	    seq[i - start] = 'A';
	    break;
	case 1:
	    seq[i - start] = 'C';
	    break;
	case 2:
	    seq[i - start] = 'G';
	    break;
	case 3:
	    seq[i - start] = 'T';
	    break;
	}
	//cerr << "  Decoding character " << i << ", '" << seq[i - start] << "'\n";
    }

    seq[i - start] = '\0';

    return seq;
}

char *DNASeq::seq(char *seq) {
    unsigned char c;
    int i;
    int index = 0;

    for (i = 0; i < this->size; i++) {
	if (i > 0 && i % bases_per_byte == 0) index++;

	//cerr << "s[" << index << "," << i % bases_per_byte << "] == " << (int)this->s[index] << "\n";

	switch (i % bases_per_byte) {
	case 0:
	    c = this->s[index] & 0xC0; // 11000000
	    c >>= 6;
	    break;
	case 1:
	    c = this->s[index] & 0x30; // 00110000
	    c >>= 4;
	    break;
	case 2:
	    c = this->s[index] & 0xC;  // 00001100
	    c >>= 2;
	    break;
	case 3:
	    c = this->s[index] & 0x3;  // 00000011
	    break;
	}

	switch (c) {
	case 0:
	    seq[i] = 'A';
	    break;
	case 1:
	    seq[i] = 'C';
	    break;
	case 2:
	    seq[i] = 'G';
	    break;
	case 3:
	    seq[i] = 'T';
	    break;
	}
	//cerr << "  Decoding character " << i << ", '" << seq[i] << "'\n";
    }

    seq[i] = '\0';

    return seq;
}

char *DNASeq::seq() {
    unsigned char c;
    int i;
    int index = 0;

    char *seq = new char[this->size + 1];

    for (i = 0; i < this->size; i++) {
	if (i > 0 && i % bases_per_byte == 0) index++;

	//cerr << "s[" << index << "," << i % bases_per_byte << "] == " << (int)this->s[index] << "\n";

	switch (i % bases_per_byte) {
	case 0:
	    c = this->s[index] & 0xC0; // 11000000
	    c >>= 6;
	    break;
	case 1:
	    c = this->s[index] & 0x30; // 00110000
	    c >>= 4;
	    break;
	case 2:
	    c = this->s[index] & 0xC;  // 00001100
	    c >>= 2;
	    break;
	case 3:
	    c = this->s[index] & 0x3;  // 00000011
	    break;
	}

	switch (c) {
	case 0:
	    seq[i] = 'A';
	    break;
	case 1:
	    seq[i] = 'C';
	    break;
	case 2:
	    seq[i] = 'G';
	    break;
	case 3:
	    seq[i] = 'T';
	    break;
	}
	//cerr << "  Decoding character " << i << ", '" << seq[i] << "'\n";
    }

    seq[i] = '\0';

    return seq;
}
