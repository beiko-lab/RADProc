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
vector<pair<string, int> > dist;
vector<int> cov;
vector<int> sample_ids;
set<string> pop_ids;
set<string> seed_matches;


Tag();
~Tag();
int  add_dist(const string seq, const int dist);
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
