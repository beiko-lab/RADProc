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

#include "tags.h"

Tag::Tag()  { 
    this->id         = 0;
   // this->cov        = 0;
    this->len        = 0;
    this->seq        = "";
    this->ascii_value = 0;
}

Tag::~Tag() { 

}
int Tag::add_dist(const string seq, const int dist) {
    //
    // Store the ID and distance as a pair, ID in the first position,
    // dist in the second.
    //
    pair<string, int> p(seq, dist);
    this->dist.push_back(p);

    return 0;
}


CTag::CTag()  { 
    this->id         = 0;
    this->len        = 0;
    //this->con        = "";
	//this->pop_count  = 0;
}

CTag::~CTag() { 
	
}

bool compare_pair_snp(pair<string, SNP *> a, pair<string, SNP *> b) {
    return (a.second->col < b.second->col);
}

int CTag::merge_snps(CTag *matched_tag)
{

 vector<SNP *>::iterator i;
    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;

    vector<pair<string, SNP *> > merged_snps;
    set<string> merged_alleles;
    set<string>::iterator s;
    SNP *csnp;

    for (i = this->snps.begin(); i != this->snps.end(); i++)
    {
	columns[(*i)->col] = make_pair("catalog", *i);
	 
	
	
	}

    for (i = matched_tag->snps.begin(); i != matched_tag->snps.end(); i++) {
	//
	// Is this column already represented from the previous sample?
	//
	if (columns.count((*i)->col)) {
	    csnp = columns[(*i)->col].second;

	    //
	    // If this is a new allele for this nucleotide, add it to the catalog SNP.
	    //
	    bool rank_1_exists = false;
	    bool rank_2_exists = false;

	    if ((*i)->rank_1 == csnp->rank_1 ||
		(*i)->rank_1 == csnp->rank_2 ||
		(*i)->rank_1 == csnp->rank_3 ||
		(*i)->rank_1 == csnp->rank_4) {
		rank_1_exists = true;
	    }
	    if ((*i)->rank_2 == csnp->rank_1 ||
		(*i)->rank_2 == csnp->rank_2 ||
		(*i)->rank_2 == csnp->rank_3 ||
		(*i)->rank_2 == csnp->rank_4) {
		rank_2_exists = true;
	    }

	    if (rank_1_exists == false) {
		if (csnp->rank_3 == 0)
		    csnp->rank_3 = (*i)->rank_1;
		else 
		    csnp->rank_4 = (*i)->rank_1;
	    }
	    if (rank_2_exists == false) {
		if (csnp->rank_3 == 0)
		    csnp->rank_3 = (*i)->rank_2;
		else 
		    csnp->rank_4 = (*i)->rank_2;
	    }

	    columns[(*i)->col] = make_pair("both", csnp);
	    
	    
	} else {
	    columns[(*i)->col] = make_pair("sample", *i);
	    
	}
    }

    for (c = columns.begin(); c != columns.end(); c++) 
	merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
     sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp);

    //
    // If the catalog tag has no defined alleles, create a matching haplotype
    // from the consensus sequence before merging in the new alleles.
    //
    string allele, new_allele;
    int    pos;

    if (this->alleles.size() == 0) {
	char c;
	new_allele = "";
	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    csnp = k->second;
	    c    = this->con[k->second->col];

	    new_allele += (csnp->col > this->len - 1) ? 'N' : c;

	    if (csnp->col > this->len - 1) continue;

	    if (c != csnp->rank_1 &&
		c != csnp->rank_2 &&
		c != csnp->rank_3 &&
		c != csnp->rank_4) {

		if (csnp->rank_3 == 0)
		    csnp->rank_3 = c;
		else 
		    csnp->rank_4 = c;
	    }
	}

	if (new_allele.length() > 0)
	    merged_alleles.insert(new_allele);
    }
  
    //
    // Merge the alleles accounting for any SNPs added from either of the two samples.
    //
    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    //
	    // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
	    // sequence to account for it in the allele string.
	    //
	    if (k->first == "sample") {
	    
		new_allele += k->second->col > this->len - 1 ? 'N' : this->con[k->second->col];
	    } else {
                new_allele += allele[pos];
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
	
	
    }

    for (j = matched_tag->alleles.begin(); j != matched_tag->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    if (k->first == "catalog") {
	     
	  
	     
		new_allele += k->second->col > matched_tag->len - 1 ? 'N' : matched_tag->con[k->second->col];
	    } else {
		new_allele += allele[pos];
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
	
	
	
    }

    //
    // If the matching tag being merged into the catalog had no called SNPs
    // create alleles from the consensus sequence and check that catalog SNP
    // objects contain all the nucleoties.
    //
    if (matched_tag->alleles.size() == 0) {
	char c;
	new_allele = "";
	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    csnp = k->second;
	    c    = matched_tag->con[k->second->col];

	    new_allele += (csnp->col > matched_tag->len - 1) ? 'N' : c;

	    if (csnp->col > matched_tag->len - 1) continue;

	    if (c != csnp->rank_1 &&
		c != csnp->rank_2 &&
		c != csnp->rank_3 &&
		c != csnp->rank_4) {

		if (csnp->rank_3 == 0)
		    csnp->rank_3 = c;
		else 
		    csnp->rank_4 = c;
	    }
	}

	if (new_allele.length() > 0)
	{
	    merged_alleles.insert(new_allele);
	    
	   
	
	    
	 }   
	    
    }

    // //
    // // If the newly merged alleles contain Ns due to different sequence lengths,
    // // check if we can reduce the alleles as one of the longer allele haplotypes
    // // may fully encompass a shorter allele haplotype that has been padded with Ns.
    // //
    // if (require_uniq_haplotypes) this->reduce_alleles(merged_alleles);

    //
    // Update the catalog entry's list of SNPs and alleles
    //
    this->snps.clear();

    for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	SNP *snp    = new SNP;
	snp->col    = (*k).second->col;
	snp->type   = (*k).second->type;
	snp->lratio = 0.0;
	snp->rank_1 = (*k).second->rank_1;
	snp->rank_2 = (*k).second->rank_2;
	snp->rank_3 = (*k).second->rank_3;
	snp->rank_4 = (*k).second->rank_4;

	this->snps.push_back(snp);

	if (k->first == "catalog" || k->first == "both")
	    delete k->second;
    }
     
     
     this->alleles.clear();
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
	this->alleles[*s] = 0;
 
    }

return 1;


}

