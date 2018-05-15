// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __GZFASTA_H__
#define __GZFASTA_H__

#ifdef HAVE_LIBZ

#include <errno.h>
#include <zlib.h>
#include "input.h"

class GzFasta: public Input {
    gzFile gz_fh;
    string buf;

 public:
    GzFasta(const char *path) : Input() { 
	this->gz_fh = gzopen(path, "rb");
	if (!this->gz_fh) {
	    cerr << "Failed to open gzipped file '" << path << "': " << strerror(errno) << ".\n";
            exit(EXIT_FAILURE);
	}
	#if ZLIB_VERNUM >= 0x1240
	gzbuffer(this->gz_fh, libz_buffer_size);
	#endif
    };
    GzFasta(string path) : Input() { 
	this->gz_fh = gzopen(path.c_str(), "rb");
	if (!this->gz_fh) {
	    cerr << "Failed to open gzipped file '" << path << "': " << strerror(errno) << ".\n";
            exit(EXIT_FAILURE);
	}
	#if ZLIB_VERNUM >= 0x1240
	gzbuffer(this->gz_fh, libz_buffer_size);
	#endif
    };
    ~GzFasta() {
	gzclose(this->gz_fh);
    };
    Seq *next_seq();
    int  next_seq(Seq &);
};

Seq *GzFasta::next_seq() {
    //
    // Check the contents of the line buffer. When we finish reading a FASTA record
    // the buffer will either contain whitespace or the header of the next FAST
    // record.
    //
    while (this->line[0] != '>' && !gzeof(this->gz_fh)) {
	gzgets(this->gz_fh, this->line, max_len);
    }

    if (gzeof(this->gz_fh)) {
	return NULL;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    //
    // Initialize the Seq structure and store the FASTA ID
    //
    Seq *s = new Seq;
    s->id = new char[len + 1];
    strcpy(s->id, this->line + 1);

    //
    // Read the sequence from the file -- keep reading lines until we reach the next
    // record or the end of file.
    //
    gzgets(this->gz_fh, this->line, max_len);

    while (this->line[0] != '>' && !gzeof(this->gz_fh)) {
	len = strlen(this->line);
	if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
	if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

	this->buf += this->line;
	gzgets(this->gz_fh, this->line, max_len);
    }

    if (gzeof(this->gz_fh)) {
	len = strlen(this->line);
	if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
	if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

	this->buf += this->line;
    }

    s->seq = new char[this->buf.length() + 1];
    strcpy(s->seq, this->buf.c_str());
    this->buf.clear();

    return s;
}

int GzFasta::next_seq(Seq &s) {
    //
    // Check the contents of the line buffer. When we finish reading a FASTA record
    // the buffer will either contain whitespace or the header of the next FAST
    // record.
    //
    while (this->line[0] != '>' && !gzeof(this->gz_fh)) {
	gzgets(this->gz_fh, this->line, max_len);
    }

    if (gzeof(this->gz_fh)) {
	return 0;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    //
    // Store the FASTA ID
    //
    strcpy(s.id, this->line + 1);

    //
    // Read the sequence from the file -- keep reading lines until we reach the next
    // record or the end of file.
    //
    gzgets(this->gz_fh, this->line, max_len);

    while (this->line[0] != '>' && !gzeof(this->gz_fh)) {
	len = strlen(this->line);
	if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
	if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

	this->buf += this->line;
	gzgets(this->gz_fh, this->line, max_len);
    }

    if (gzeof(this->gz_fh)) {
	len = strlen(this->line);
	if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
	if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

	this->buf += this->line;
    }

    strcpy(s.seq, this->buf.c_str());
    this->buf.clear();

    return 1;
}

#else  // If HAVE_LIBZ is undefined and zlib library is not present.

#include "input.h"

class GzFasta: public Input {
 public:
    GzFasta(const char *path) : Input() { cerr << "Gzip support was not enabled when Stacks was compiled.\n"; };
    ~GzFasta() {};
    Seq *next_seq()      { return NULL; };
    int  next_seq(Seq &) { return 0; };
};

#endif // HAVE_LIBZ

#endif // __GZFASTA_H__
