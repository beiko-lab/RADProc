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

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

//
// Pull in the configuration variables from the configure script
//
#if HAVE_CONFIG_H
#include "config.h"
#endif

//
//
//
const unsigned int fieldw = 4;

//
// Maximum line length for parsing input files.
//
const int max_len = 1024;

//
// Maximum length of idetifiers, such as sequence IDs and chromosome names.
//
const int id_len = 255;

//
// Size to use for internal buffer size for gzipped files being read with libz.
//
const int libz_buffer_size = 1048576;

//
// Supported file types
//
enum file_type {unknown, 
		sql,     gzsql, 
		fasta,   gzfasta, 
		fastq,   gzfastq, 
		bowtie,  sam, bam, tsv, 
		bustard, phase, fastphase, beagle};

#endif
