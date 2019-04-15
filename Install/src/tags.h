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

#ifndef __TAGS_H__
#define __TAGS_H__
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;
#include <utility>
using std::pair;
using std::make_pair;
#include<iostream>
using std::cerr;
#include <algorithm>
using std::cout;

class Tag{
public:
string seq;
int id;
//int cov;
int len;
int ascii_value;
vector<pair<int, int> > dist;
vector<int> cov;
vector<int> sample_ids;
set<string> pop_ids;
set<int> seed_matches;


Tag();
~Tag();
int  add_dist(const int id, const int dist);
};

enum snp_type    {snp_type_het, snp_type_hom, snp_type_unk};
class SNP {
 public:
    snp_type type;   // Heterozygous or homozygous
    int     col;
    float    lratio;
    char     rank_1;
    char     rank_2;
    char     rank_3;
    char     rank_4;

    SNP() {
	col    = 0;
	lratio = 0.0;
	rank_1 = 0;
	rank_2 = 0;
	rank_3 = 0;
	rank_4 = 0;
    }
};


class CTag{
public:
	int id;
	string con;
	vector<SNP *> snps;
	map<string,int> alleles;
	set<pair<int, int> > sources;
	int len;
	std::vector <std::string> haps;
	map<string, set<int> > sam_per;
	map< string, std::pair <std::string, int > > allele_strings; 
	CTag();
	~CTag();
    int merge_snps(CTag *);
  
};

class MTag: public CTag{
public:
  int cat_id;   

};




#endif
