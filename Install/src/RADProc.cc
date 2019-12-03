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

//
// RADProc -- De novo locus assembler and catalog builder for RADSeq data
//


#include "RADProc.h"

string    in_file;
string    in_file_type;
string    out_path;
int       num_threads       = 1;
int       max_nwk_dist     = 4;
int       max_cat_dist     = 1;
bool      call_sec_hapl     = true;
int       min_merge_cov     = 3;
int       max_rem_dist      = -1;
double    min_sam_per     = 0.1;
int       max_stacks = 3;
int    min_depth     = 2;
int psweep =0;
double alpha              = 0.05;
double heterozygote_limit = -3.84;
double homozygote_limit   =  3.84;
int num_het_locus = 0;
int num_hom_locus = 0;
int num_2_alleles = 0;
int num_3_alleles = 0;
int num_more_3_alleles =0;

int tot_loci_more_3_snps = 0; 
int id = 0;
std::ofstream log_file;



void help() {
    std::cerr << "RADProc" << "1.0" << "\n"
              << "RADProc -t infile_type -f file_path [-o path][-a psweep][-M max_dist] [-m min_cov][-n max_cat_dist][-p num_threads] [-x max_stacks][-S min_sam][-D min_depth][-h]" << "\n"
              << "  t: Input file Type. Supported types: fasta, fastq, gzfasta, or gzfastq.\n"
              << "  f: Input file path.\n"
	          << "  o: Output path to write results.\n"
	          << "  a: Enable parameter sweep mode.\n"
	          << "  M: Maximum distance (in nucleotides) allowed between stacks to form network.\n"
	          << "  m: Minimum coverage depth.\n"
	          << "  n: Maximum distance (in nucleotides) allowed between catalog loci to merge.\n"
              << "  p: Enable parallel execution with num_threads threads.\n"
              << "  x: Maximum number stacks per locus.\n"
              << "  S: Minimum sample percentage.\n"
              << "  D: Minimum average coverage depth.\n"
	          << "  h: Display this help messsage.\n";


    exit(0);
}

int parse_command_line(int argc, char* argv[]) {
    int c;
    int count;
    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
	    {"file",         required_argument, NULL, 'f'},
        {"infile_type",  required_argument, NULL, 't'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"max_dist",     required_argument, NULL, 'M'},
	    {"min_cov",     required_argument, NULL, 'm'},
	    {"max_cat_dist",     required_argument, NULL, 'n'},
	    {"num_threads",  required_argument, NULL, 'p'},
	    {"max_stacks",  required_argument, NULL, 'x'},
	    {"min_sam",  required_argument, NULL, 'S'},
	    {"min_depth",  required_argument, NULL, 'D'},
	    {"Parameter search",  no_argument, NULL, 'a'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
        count++;
	c = getopt_long(argc, argv, "haf:t:o:p:x:M:m:n:S:D:", long_options, &option_index);
	// Detect the end of the options.
	
        if (c == -1)
           {
            if ( count == 1)
             {
               help();
	       abort();  
             }
            else
	    break;
           }
       
	switch (c) {
	case 'h':
	    help();
	    break;
        case 't':
            if (strcmp(optarg, "fasta") == 0)
                in_file_type = "fasta";
            else if (strcmp(optarg, "fastq") == 0)
                in_file_type = "fastq";
	    else if (strcasecmp(optarg, "gzfasta") == 0)
                in_file_type = "gzfasta";
	    else if (strcasecmp(optarg, "gzfastq") == 0)
                in_file_type = "gzfastq";
            else
                in_file_type = "unknown";
     	case 'f':
	    in_file = optarg;
	    break;
	case 'o':
	    out_path = optarg;
	    break;
	case 'M':
	    max_nwk_dist = atoi(optarg);
	    break;
	case 'm':
	    min_merge_cov = atoi(optarg);
	    break; 
	case 'n':
	    max_cat_dist = atoi(optarg);
	    break; 
	case 'S':
	    min_sam_per = atof(optarg);
	case 'D':
	    min_depth = atof(optarg);    
	    break;       
	case 'p':
	    num_threads = atoi(optarg);
	    break;
	case 'x':
	      max_stacks = atoi(optarg);
	      break;   
	case 'a':
	    psweep = 1; //atoi(optarg);
	    break;    
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
	default:
	    cerr << "Unknown command line option '" << (char) c << "'\n";
	    help();
	    abort();
	}
    }

    if (in_file.length() == 0 || in_file_type == "unknown") {
	cerr << "You must specify an input file of a supported type.\n";
	help();
    }
    
    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";
    return 0;
}

int add_sum(string allele) {

          int ascii_val=0, val=0;
          for(int ascii=0; ascii < allele.length(); ascii++){
	       
	         if ( allele.at(ascii) == 'A') val=1;
             else if ( allele.at(ascii) == 'T') val=2;
             else if ( allele.at(ascii) == 'G') val=3;
             else  val=4;
             //tag->ascii_value = tag->ascii_value + int(tag -> seq.at(j));
             ascii_val= ascii_val + val;
	         // ascii_val = ascii_val+int(alelle_s.at(ascii));
	       }
	       //pair<string,int> p(allele, ascii_val);
           //this->alleles.push_back(p);
    return  ascii_val;          

}	

bool compare_pair(pair<char, int> a, pair<char, int> b) {
    return (a.second > b.second);
}


snp_type
call_multinomial_snp(vector<SNP *> &ctag_snps, vector<SNP *> &snps, int col, map<char, int> &n, bool record_snps) 
{   
    vector<pair<char, int> > nuc;
    map<char, int>::iterator i;
    
    int total = 0;
    for (i = n.begin(); i != n.end(); i++) {
	if (i->first != 'N') {
	    total += i->second;
	    nuc.push_back(make_pair(i->first, i->second));
	}
    }

    sort(nuc.begin(), nuc.end(), compare_pair);

    //
    // If this column was simply uncalled Ns, return.
    //
    if (nuc[0].second == 0) {
	if (record_snps) {
	    SNP *snp = new SNP;
	    snp->type   = snp_type_unk;
	    snp->col    = col;
	    snp->lratio = 0;
	    snp->rank_1 = 'N';
	    snp->rank_2 = '-';
        //delete snp;
	    snps.push_back(snp);
	}
	return snp_type_unk;
    }

    double nuc_1   = nuc[0].second;
    double nuc_2   = nuc[1].second;
    double nuc_3   = nuc[2].second;
    double nuc_4   = nuc[3].second;
    double l_ratio = 0;

    l_ratio = (nuc_1 * log(nuc_1 / total));

    if (total - nuc_1 > 0)
    	l_ratio += ((total - nuc_1) * log((total - nuc_1) / (3 * total)));

    if (nuc_1 + nuc_2 > 0)
	l_ratio -= ((nuc_1 + nuc_2) * log((nuc_1 + nuc_2) / (2 * total)));

    if (nuc_3 + nuc_4 > 0)
	l_ratio -= ((nuc_3 + nuc_4) * log((nuc_3 + nuc_4) / (2 * total)));

    l_ratio *= 2;

    snp_type res;

    if (l_ratio <= heterozygote_limit) {
	//
        // This locus is a heterozygote.
	//
	if (record_snps) {
	    SNP *snp = new SNP;
	    snp->type   = snp_type_het;
	    snp->col    = col;
	    snp->lratio = l_ratio;
	    snp->rank_1 = nuc[0].first;
	    snp->rank_2 = nuc[1].first;
       
           // cout<<snp->type<<"\n";

	   snps.push_back(snp);
	   //delete snp;
	   SNP *ctag_snp = new SNP;
	    ctag_snp->type   = snp_type_het;
	    ctag_snp->col    = col;
	    ctag_snp->lratio = l_ratio;
	    ctag_snp->rank_1 = nuc[0].first;
	    ctag_snp->rank_2 = nuc[1].first;

	   ctag_snps.push_back(ctag_snp);
	   //delete ctag_snp;
	}
	res = snp_type_het;

    } else if (l_ratio >= homozygote_limit) {
	//
        // This locus is a homozygote.
	//
	if (record_snps) {
	    SNP *snp = new SNP;
	    snp->type   = snp_type_hom;
	    snp->col    = col;
	    snp->lratio = l_ratio;
	    snp->rank_1 = nuc[0].first;
	    snp->rank_2 = '-';
 
	    snps.push_back(snp);
	    //delete snp;
	
	    
	}
	res = snp_type_hom;

    } else {
	//
        // Unknown whether this is a heterozygote or homozygote.
	//
	if (record_snps) {
	    SNP *snp = new SNP;
	    snp->type   = snp_type_unk;
	    snp->col    = col;
	    snp->lratio = l_ratio;
	    snp->rank_1 = nuc[0].first;
	    snp->rank_2 = nuc[1].second > 0 ? nuc[1].first : '-';

	    snps.push_back(snp);
	    //delete snp;
	    

	}

	res = snp_type_unk;
    }
    
    

    return res;
}	

bool compare_pair_snp1(pair<string, SNP *> a, pair<string, SNP *> b) {
    return (a.second->col < b.second->col);
}

int merge_allele(CTag *locus, SNP *snp) {
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;
    vector<SNP *>::iterator i;
    SNP *lsnp;

    for (i = locus->snps.begin(); i != locus->snps.end(); i++)
    {
     //if ((*i)->type == snp_type_het)
	columns[(*i)->col] = make_pair("sample", *i);
	
	
	}

    if (columns.count(snp->col)) {
	lsnp = columns[snp->col].second;

	//
	// If this is a new allele for this nucleotide, add it to the catalog SNP.
	//
	bool rank_1_exists = false;
	bool rank_2_exists = false;

	if (snp->rank_1 == lsnp->rank_1 ||
	    snp->rank_1 == lsnp->rank_2 ||
	    snp->rank_1 == lsnp->rank_3 ||
	    snp->rank_1 == lsnp->rank_4) {
	    rank_1_exists = true;
	}
	if (snp->rank_2 == lsnp->rank_1 ||
	    snp->rank_2 == lsnp->rank_2 ||
	    snp->rank_2 == lsnp->rank_3 ||
	    snp->rank_2 == lsnp->rank_4) {
	    rank_2_exists = true;
	}

	if (rank_1_exists == false) {
	    if (lsnp->rank_3 == 0)
		lsnp->rank_3 = snp->rank_1;
	    else 
		lsnp->rank_4 = snp->rank_1;
	}
	if (rank_2_exists == false) {
	    if (lsnp->rank_3 == 0)
		lsnp->rank_3 = snp->rank_2;
	    else 
		lsnp->rank_4 = snp->rank_2;
	}

	columns[snp->col] = make_pair("both", lsnp);
	
    } else {
	columns[snp->col] = make_pair("merge", snp);
	
	
    }

    vector<pair<string, SNP *> > merged_snps;

    for (c = columns.begin(); c != columns.end(); c++) 
	merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp1);

    //
    // Modify any existing alleles to account for this new SNP. If there are not any alleles, 
    // create new ones.
    //
    stringstream sallele;
    set<string> merged_alleles;
    string allele, new_allele;
    int pos;

    if (locus->alleles.size() == 0) {
        sallele << locus->con[snp->col];
        merged_alleles.insert(sallele.str());
    }

    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;    

    for (j = locus->alleles.begin(); j != locus->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

        // cerr << "Allele length: " << allele.size() << "\n";

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    //
	    // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
	    // sequence to account for it in the allele string.
	    //
	    if ((*k).first == "merge") {
		new_allele += locus->con[(*k).second->col];
                // cerr << "  Adding char '" << locus->con[k->second->col] << "' from consensus position " << (*k).second->col << "\n";
	    } else {
		new_allele += allele[pos];
                // cerr << "  Adding char '" << allele[pos] << "' from allele position " << pos << "\n";
		pos++;
	    }
	    
	  
	}

	merged_alleles.insert(new_allele);
    }

    set<string>::iterator s;

    locus->alleles.clear();
    string s1 = "";
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
     s1= *s;
     if (s1.find("N") == std::string::npos) {
     //cout << s1 
	 locus->alleles[s1] = 0;
	 
	 }
    }

    return 1;
}

int
characterize_mismatch_snps(CTag *catalog_tag, CTag *query_tag)
{
    set<int> snp_cols;
    uint i;
    for (i = 0; i < catalog_tag->snps.size(); i++)
    {
    //if (catalog_tag->snps[i]->type == snp_type_het)
	snp_cols.insert(catalog_tag->snps[i]->col);
	
	}
    for (i = 0; i < query_tag->snps.size(); i++)
    {
    //if (query_tag->snps[i]->type == snp_type_het)
	snp_cols.insert(query_tag->snps[i]->col);
	
	}
	

    //
    // For each mismatch found, create a SNP object
    //
    const char *c        = catalog_tag->con.c_str();
    const char *c_beg    = c;
    const char *c_end    = c + strlen(c);
    const char *q        = query_tag->con.c_str();
    const char *q_beg    = q;
    const char *q_end    = q + strlen(q);

    i = 0;
    while (c < c_end && q < q_end) {

	if (snp_cols.count(i) == 0 && (*c != *q) && (*c != 'N' && *q != 'N')) {

            // cerr << "Adding a new SNP at position " << c - c_beg << ", " << *c << "/" << *q << "\n";
            SNP *s = new SNP;
	    s->type   = snp_type_het;
            s->col    = c - c_beg;
            s->lratio = 0;
            s->rank_1 = *c;
            s->rank_2 = *q;

            merge_allele(catalog_tag, s);
            merge_allele(query_tag, s);
            
            catalog_tag->snps.push_back(s);

            s = new SNP;
	        s->type   = snp_type_het;
            s->col    = q - q_beg;
            s->lratio = 0;
            s->rank_1 = *q;
            s->rank_2 = *c;
            
            query_tag->snps.push_back(s);
        }
	c++;
	q++;
	i++;
    }
 
    return 1;
}

int call_primary_consensus( map<string, vector<string> > &merged, map<string, Tag *> ptags, map<string,string> &consensus, int sample_id) {

std::map<string, std::vector<string> >::iterator it;
vector<string> keys;
 vector<int>::iterator vec_it;
 for (it = merged.begin(); it != merged.end(); it++) 
       {
       // if ( (it ->second.size() >= 1) && (it ->second.size() <= 3) )
        if (it ->second.size() <= 4)
        {
    	keys.push_back(it->first);
    	
    	}
        
       }
cerr << "\n" << "Calling consensus for " << sample_id << "\n\n";
	std::vector<string> merged_tag;
	std::vector<string>::iterator merged_tag_it;
	map<int, map <int, vector<string> > > reads_map;
	map<int, vector<string> >::iterator reads_map_it;
	Tag       *ptag;
	cerr << "\n" << "Number of Loci " << keys.size() << "\n\n";
	
	
for (int i = 0; i < (int)keys.size(); i++) {
    	
          //if (i % 1000 == 0) cerr << "Processing  Merged-Tag " << i << "       \r";
        
    	    merged_tag = merged[keys[i]];
    	   // cout <<  keys[i] << "\n";

    	    vector<string>::iterator j;
    	    vector<string>  reads;
    	     vector<string>  temp_reads;
            map<string,int> alleles;
            char    base;
            vector<SNP *>::iterator snp;
            string allele;
            int height = 0;
            int length = 0;

			
		    int l = 0;
    	    for (j = merged_tag.begin(); j != merged_tag.end(); j++) {    
    		     ptag = ptags[*j]; 
                length = ptag->len;
                int pos=0;
                for( vec_it = ptag->sample_ids.begin(); vec_it != ptag->sample_ids.end(); vec_it++)
                {
                   if (sample_id == *vec_it) break;
                   pos++;
                }
            
                
                if ( ptag->cov.size() == pos) continue;
                
                for ( int k =0; k < ptag->cov.at(pos);k++)
                { 
                 reads.push_back(ptag->seq);
                 temp_reads.push_back(ptag->seq);  
                } 
              reads_map[i][l] =temp_reads;
              temp_reads.clear();
              
              l++;  
                
    	    }
    	    
            height = reads.size();
    	    int row, col;
    	    string con;
    	    map<char, int> nuc;
    	    map<char, int>::iterator max, n;
            string d;

    	    for (col = 0; col < length; col++) {
    		nuc['A'] = 0; 
    		nuc['G'] = 0;
    		nuc['C'] = 0;
    		nuc['T'] = 0;

    		for (row = 0; row < height; row++) {
    		    d = reads[row];
		    if (nuc.count(d[col]))
			nuc[d[col]]++;
    		}
    		//
    		// Find the base with a plurality of occurances and call it.
    		//
    		max = nuc.end();

    		for (n = nuc.begin(); n != nuc.end(); n++) {

    		    if (max == nuc.end() || n->second > max->second)
    			max = n;
    		}
    		con += max->first;

    	    }
    	    
    	consensus[keys[i]]=con;    
    
   }
return 0;
}



int call_consensus(map<int, vector<int> > &merged, map<int, Tag *> ptags,map<int, CTag *> &ctags,map<int,set<int> > &ctag_alleles_map, int sample_id, string sample_name, string path, int cov, string param,  ofstream &stats_fh) 
{
  
    std::ofstream snp_s,al,mods, tags, cat_snps, cat_all, mat;
    std::map<int, std::vector<int> >::iterator it;
    vector<int> keys;
    map<int, map<string,int> > alleles_map;   
    map<int, map<string,int> >::iterator als_it;
    map<string,int>::iterator al_it;
    map<int, vector<SNP *> > snps_map;
    map<int, vector<SNP *> >::iterator snps_it;
    map<int, vector<int> >::iterator merged_it;
    
	std::map<int, std::map< string,int > >::iterator ctags_pop_count_it;
    vector<SNP*>::iterator snps_it1;
    std::vector<std::pair <std::string, int > >::iterator cat_alleles_it;
  
    
  
  
    string a;
    string snp_file;
    string all_file;
    string mod_file;
    string tags_file;
    string matches_file;
    string cat_snps_file;
    string cat_all_file;
    
    
    

    
	 map<int,int> num_snps;
    map<int,int>::iterator num_snps_it;
	
	vector<SNP *> ctag_snps;
	int val;
    int ascii_val = 0;
    std::map<string, std::vector<std::pair <std::string, int > > > mtag_alleles;	 

   std::map<int, set<int> >::iterator matches_it;
   set<int>::iterator matches_it1;
   
    
    snp_file = path +sample_name+".snps.tsv";
    all_file = path +sample_name+".alleles.tsv";
    mod_file = path +sample_name+".models.tsv";
    tags_file = path +sample_name+".tags.tsv";
    //matches_file = path +sample_name+".matches.tsv";
     
    
    int tot_num_snps=0;
    int tot_poly_loci=0;
    int poly_flag=0;
    
    vector<int>::iterator vec_it;
    map<int,vector<int> > ctag_matches;
    map<int,set<int> >::iterator ctag_alleles_map_it;
    

    
     vector<pair<int, int> >::iterator dist_it; 
     
   
     
     
     
     
   for (it = merged.begin(); it != merged.end(); it++) 
       {
        if (it->second.size() > 0) 
    
        {
    	keys.push_back(it->first);
    	
    	}
        
       }
    stringstream log;   
     time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);   
    snp_s.open(snp_file.c_str());
    al.open(all_file.c_str());
    mods.open(mod_file.c_str());
    tags.open(tags_file.c_str());
   // mat.open(matches_file.c_str());
    

    log << "# RAD2POP version 1.0 generated on " << date << "\n";
    
   
 
        tags << log.str();
        mods << log.str();
        snp_s << log.str();
        al << log.str();
    int i=0, id;
    int count=0;
     cerr << "\n" << "Calling consensus for " << sample_name << "\n\n";
	std::vector<int> merged_tag;
	std::vector<int>::iterator merged_tag_it;
	map<int, map <int, vector<string> > > reads_map;
	map<int, vector<string> >::iterator reads_map_it;
	Tag       *ptag;
    int al_count = 0, snp_count=0, stack_count=0, loci_count =0;
        
    int num_loci =0, empty_merge = 0;
    	for (i = 0; i < (int)keys.size(); i++) {
    	
    	   stack_count=0;
    	
           if (i % 1000 == 0) cerr << "Processing  Merged-Tag " << i << "       \r";
            vector<SNP *> snps; 
            set<pair<int,string> > alleles_str;
            set<pair<int,string> >::iterator alleles_str_it;
            set<int> con_str;
            set<int>::iterator con_str_it;
    	    merged_tag = merged[keys[i]];
    	    vector<int>::iterator j;
    	    vector<string>  reads;
    	     vector<string>  temp_reads;
            map<string,int> alleles;
            map<string, int> ptag_map;
            char    base;
            vector<SNP *>::iterator snp;
            string allele;
            int height = 0;
            int length = 0;
            CTag *ctag;
	        CTag *ctag1;
			ctag = new CTag();
		    int l = 0;
		     for (j = merged_tag.begin(); j != merged_tag.end(); j++) {     
    		     ptag = ptags[*j];
                 ptag_map[ptag->seq]=ptag->id;
                length = ptag->len;
                int pos=0;
                for( vec_it = ptag->sample_ids.begin(); vec_it != ptag->sample_ids.end(); vec_it++)
                {
                   if (sample_id == *vec_it) break;
                   pos++;
                }
            
                
                if ( ptag->cov.size() == pos) continue;
                
                if ( ptag->cov.at(pos) >= cov)
                {
                stack_count++;

                }
                
                for ( int k =0; k < ptag->cov.at(pos);k++)
                { 
                 reads.push_back(ptag->seq);
                 temp_reads.push_back(ptag->seq);  
                } 
              reads_map[i][l] =temp_reads;
              temp_reads.clear();
              
              l++;  
                
    	    }
    	    
            height = reads.size();
    	    int row, col;
    	    string con;
    	    map<char, int> nuc;
    	    map<char, int>::iterator max, n;
            string d;

    	    for (col = 0; col < length; col++) {
    		nuc['A'] = 0; 
    		nuc['G'] = 0;
    		nuc['C'] = 0;
    		nuc['T'] = 0;

    		for (row = 0; row < height; row++) {
    		    d = reads[row];
		    if (nuc.count(d[col]))
			nuc[d[col]]++;
    		}
    		//
    		// Find the base with a plurality of occurances and call it.
    		//
    		max = nuc.end();

    		for (n = nuc.begin(); n != nuc.end(); n++) {

    		    if (max == nuc.end() || n->second > max->second)
    			max = n;
    		}
    		con += max->first;
            
		  call_multinomial_snp(ctag->snps,snps, col, nuc, true);

    	    }
  
           snps_map[i]= snps;
            int ascii_val=0;
           
        uint snp_cnt;
           for (row = 0; row < height; row++) {
	   allele.clear();

	   snp_cnt = 0;

	   for (snp = snps.begin(); snp != snps.end(); snp++) {
	   if ((*snp)->type != snp_type_het) continue;

	   snp_cnt++;
          
          
           d  = reads[row];
	       base = d[(*snp)->col];

	    if (base == (*snp)->rank_1 || base == (*snp)->rank_2) 
		allele += base;
	    else
		break;
	    }
  
	   if (snp_cnt > 0 && allele.length() == snp_cnt)
	   {
	    alleles[allele]++;
	    
	   ctag->alleles[allele]++;
	    pair<int, string> p( ptag_map[d], allele);
	    alleles_str.insert(p);
	    }
	    
	    if (snp_cnt == 0)
           {   
               d = reads[0];
               con_str.insert(ptag_map[d]);
        
           }
            
	    
           }
           
ptag_map.clear();
         
         for(al_it = alleles.begin(); al_it != alleles.end(); al_it++)
         {
         
         al << "0" << "\t" << sample_id-1 << "\t" << i <<"\t"<<al_it->first <<"\t"<< ((double)al_it->second/height)*100 <<"\t" <<al_it->second <<"\n";
	
         }
         
         alleles_map[i]=alleles;
         
         if ( alleles.size() > 0) 
         {
         snp_count = snp_count + snp_cnt;
         al_count++;
         }
     if ( sample_id == 1)
          {    
          ctag->id = ctags.size();
          ctag->con = con;
          ctag->len = con.length();
          for ( con_str_it = con_str.begin(); con_str_it != con_str.end(); con_str_it++)
             {
              ctag_alleles_map[*con_str_it].insert(ctag->id); 
                     
              }
         
          for ( alleles_str_it = alleles_str.begin(); alleles_str_it != alleles_str.end(); alleles_str_it++)
             {
             ctag_alleles_map[alleles_str_it->first].insert(ctag->id);         
             }
             ctag->sources.insert(make_pair(sample_id-1, i));
             ctag->sam_per[sample_name.substr(0,3)].insert(sample_id-1);
           ctags[ctags.size()] = ctag; 
         
          
         }
        
        if ( sample_id > 1)

          {
             std::map<int, set<int> > matches;
  
             for ( con_str_it = con_str.begin(); con_str_it != con_str.end(); con_str_it++)
             {
            
             ptag = ptags[*con_str_it];
               
             for( dist_it = ptag->dist.begin(); dist_it != ptag->dist.end(); dist_it++)
                {    
                
                   if (dist_it->second <= max_cat_dist)
                   {  
                   
                    ctag_alleles_map_it = ctag_alleles_map.find(dist_it->first);
                   
                    if (ctag_alleles_map_it != ctag_alleles_map.end())
                    {           
                     matches[dist_it->second].insert( *ctag_alleles_map[dist_it->first].begin());
                    
                    }
                   
                   }
                   
                   
                }
             
            
            }

         for ( alleles_str_it = alleles_str.begin(); alleles_str_it != alleles_str.end(); alleles_str_it++)
             {
             ptag = ptags[alleles_str_it->first];
          
             for( dist_it = ptag->dist.begin(); dist_it != ptag->dist.end(); dist_it++)
                {              
                    if (dist_it->second <= max_cat_dist)
                   {
                    ctag_alleles_map_it = ctag_alleles_map.find(dist_it->first);
                    if (ctag_alleles_map_it != ctag_alleles_map.end())
                    {
                   matches[dist_it->second].insert( *ctag_alleles_map[dist_it->first].begin());

             
                   }
                   
                   }
                   
                  
                }
             
             }
            
           if (( matches.begin()->second.size() == 0) || ( matches.begin()->second.size() > 1000000) )
           {
                                               
             ctag->id = ctags.size();
             ctag->con = con;
             ctag->len = con.length();
             ctag->sources.insert(make_pair(sample_id-1, i));
             ctag->sam_per[sample_name.substr(0,3)].insert(sample_id-1);
             for( merged_tag_it = merged_tag.begin(); merged_tag_it != merged_tag.end(); merged_tag_it++)
            {
             ptag = ptags[*merged_tag_it];
             vec_it = find(ptag->sample_ids.begin(), ptag->sample_ids.end(), sample_id);
             if ( vec_it  != ptag->sample_ids.end())
             ctag_alleles_map[*merged_tag_it].insert(ctag->id);
                                  
            }        
           ctags[ctags.size()] = ctag; 
           } 
           
        
             
           else if ( matches.begin()->second.size() == 1)
           {
            for( merged_tag_it = merged_tag.begin(); merged_tag_it != merged_tag.end(); merged_tag_it++)
            {
            ctag_alleles_map_it = ctag_alleles_map.find(*merged_tag_it);
              if (ctag_alleles_map_it == ctag_alleles_map.end())
              {
              ptag = ptags[*merged_tag_it];
             vec_it = find(ptag->sample_ids.begin(), ptag->sample_ids.end(), sample_id);
             if ( vec_it  != ptag->sample_ids.end())
            ctag_alleles_map[*merged_tag_it].insert(*matches.begin()->second.begin());  
              }          
                                  
            }
            ctag->con = con;
            ctag->len = con.length();
            ctag1 = ctags[*matches.begin()->second.begin()];
            characterize_mismatch_snps(ctag1,ctag);
            ctag1->merge_snps(ctag);
           ctag1->sources.insert(make_pair(sample_id-1, i));
           ctag1->sam_per[sample_name.substr(0,3)].insert(sample_id-1);
            
           }  
   
            matches.clear();
                  
         }     
      
   }
   

  for (snps_it = snps_map.begin(); snps_it != snps_map.end(); snps_it++)
	{
         int snp_cnt = 0;
        
         poly_flag=0;
        
      if ( snps_it->second.size() > 0)
      {
       num_loci++;
     mods<< "0" << "\t" << sample_id-1<< "\t" << snps_it->first<<"\t"  << ""<< "\t"  << 0 << "\t"  << "+" << "\t" << "consensus\t" << "\t\t";
    tags<< "0" << "\t" << sample_id-1<< "\t" << snps_it->first<<"\t"  << ""<< "\t"  << 0 << "\t"  << "+" << "\t" << "consensus\t" << "\t\t";

    for (snps_it1 = snps_it->second.begin();snps_it1 != snps_it->second.end();snps_it1++)
	{
       mods << (*snps_it1)->rank_1 ;
       tags << (*snps_it1)->rank_1 ;
    }  
    
     mods << "\t" << "0" << "\t" << "0" << "\t" << "0" <<"\t" << "0"<<"\n"; 
     tags << "\t" << "0" << "\t" << "0" << "\t" << "0" <<"\t" << "0"<<"\n";
     mods<< "0" << "\t" << sample_id-1<< "\t" << snps_it->first<<"\t"  << ""<< "\t"  << 0 << "\t"  << "+" << "\t" << "model\t"    << "\t\t" ;
     tags<< "0" << "\t" << sample_id-1<< "\t" << snps_it->first<<"\t"  << ""<< "\t"  << 0 << "\t"  << "+" << "\t" << "model\t"    << "\t\t" ;
    
    for (snps_it1 = snps_it->second.begin();snps_it1 != snps_it->second.end();snps_it1++)
	{
           
           if( (*snps_it1)->type == 0) 
         {
           tot_num_snps++;
           snp_cnt++;
           if ( snp_cnt > 3) tot_loci_more_3_snps++;
           poly_flag=1;
         }
      snp_s<< "0" << "\t" << sample_id-1<< "\t" << snps_it->first<<"\t" << (*snps_it1)->col  << "\t";  
             switch((*snps_it1)->type) {
            case 0:
                snp_s << "E\t";
                mods << "E";
                tags << "E";
                break;
            case 1:
                snp_s << "O\t";
                mods << "O";
                tags << "O";
                break;
            default:
                snp_s << "U\t";
                mods << "U";
                tags << "U";
                break;
            }  
	 snp_s << (*snps_it1)->lratio << "\t" << (*snps_it1)->rank_1  << "\t" << (*snps_it1)->rank_2 << "\n";
	delete *snps_it1;
	}
	mods << "\t\t\t\t"<<"\n"; 
	tags << "\t\t\t\t"<<"\n"; 
	num_snps[snp_cnt]++;    
    id = 0;
   
   
   for ( reads_map_it = reads_map[snps_it->first].begin(); reads_map_it != reads_map[snps_it->first].end(); reads_map_it++)
    {
              if ( reads_map_it->second.size() >= cov )
                {
                for ( int k =0; k < reads_map_it->second.size();k++)
                {
                              
                tags << "0"       << "\t"
                     << sample_id-1    << "\t"
                     << snps_it->first << "\t"
                     << "\t" // chr
                     << "\t" // bp
                     << "\t" // strand
                     << "primary\t" 
                     << id << "\t" 
                     << k << "\t" 
                     << *reads_map_it->second.begin()
                     << "\t\t\t\t\n";
                }   
                }
                else
                {
                for ( int k =0; k < reads_map_it->second.size();k++)
                {
                
                tags << "0"       << "\t"
                     << sample_id-1    << "\t"
                     << snps_it->first << "\t"
                     << "\t" // chr
                     << "\t" // bp
                     << "\t" // strand
                     << "secondary\t" 
                     << "\t" 
                     << k << "\t" 
                     << *reads_map_it->second.begin()
                     << "\t\t\t\t\n";
                }   
                id--;
                
                }
                
                 
                id++;
       }       
    
	
         if (snp_cnt > 0)  num_het_locus++;
         else num_hom_locus++;
     }   
     
     if (poly_flag == 1) tot_poly_loci++;
	}

stats_fh << sample_name << "\t" << param << "\t" << num_loci << "\t" << tot_poly_loci << "\t" << tot_num_snps << "\n"; 

cerr << "\n" << "Number of loci: " << num_loci << "\n";
cerr << "\n" << "Polymorphic loci: " << tot_poly_loci << "\n";
cerr << "\n" << "SNPs " <<  tot_num_snps  <<"\n";


return 0;
}

int dist(Tag *tag_1,Tag *tag_2, int distance) {
    int   dist  = 0;
    const char *p     = tag_1->seq.c_str();
    const char *q     = tag_2->seq.c_str();
    const char *p_end = p + tag_1->len;
    const char *q_end = q + tag_1->len;
    int len = 0;  
    //
    // Count the number of characters that are different
    // between the two sequences.
    //
    while (p < p_end && q < q_end) {
    len++;
	dist += (*p == *q) ? 0 : 1;
	p++; 
	q++;
      
	if (dist > distance)
		return -1;
	//if ( len >= tag_1->len/2) break;	
        
    }

    return dist;
}

int dist(std::string tag_1,std::string tag_2, int distance) {
    int   dist  = 0;
    const char *p     = tag_1.c_str();
    const char *q     = tag_2.c_str();
    const char *p_end = p + tag_1.length();
    const char *q_end = q + tag_1.length();
      
    //
    // Count the number of characters that are different
    // between the two sequences.
    //
    while (p < p_end && q < q_end) {
	dist += (*p == *q) ? 0 : 1;
	p++; 
	q++;
      
	if (dist > distance)
		return -1;
        
    }

    return dist;
}

int determine_kmer_length(int read_len, int dist) {
    int kmer_len, span, min_matches;

    //
    // If distance allowed between sequences is 0, then k-mer length equals read length.
    //
    if (dist == 0) 
        return read_len;

    //
    // Longer k-mer lengths will provide a smaller hash, with better key placement.
    // Increase the kmer_len until we start to miss hits at the given distance. Then
    // back the kmer_len off one unit to get the final value.
    //
    for (kmer_len = 5; kmer_len < read_len; kmer_len += 2) {
        span = (kmer_len * (dist + 1)) - 1;

        min_matches = read_len - span;

        if (min_matches <= 0) break;
    }

    if (kmer_len >= read_len) {
        cerr << "Unable to find a suitable k-mer length for matching.\n";
        exit(1);
    }

    kmer_len -= 2;

    return kmer_len;
}

int calc_min_kmer_matches(int kmer_len, int dist, int read_len, bool exit_err) {
    int span, min_matches;

    span = (kmer_len * (dist + 1)) - 1;

    min_matches = read_len - span;

    if (min_matches <= 0) {
        cerr << 
            "Warning: combination of k-mer length (" << kmer_len << ") and edit distance (" << dist << ") allows for " <<
            "sequences to be missed by the matching algorithm.\n";
    }

    if (min_matches <= 0 && exit_err)
        exit(1);
    else if (min_matches <= 0)
        min_matches = 1;

    return min_matches;
}
int
initialize_kmers(int kmer_len, int num_kmers, vector<char *> &kmers)
{
    char *kmer;

    for (int i = 0; i < num_kmers; i++) {
        kmer = new char[kmer_len + 1];
        kmers.push_back(kmer);
    }

    return 0;
}

int
generate_kmers_lazily(const char *seq, uint kmer_len, uint num_kmers, vector<char *> &kmers)
{
    char *kmer;
    const char *k = seq;

    if (num_kmers > kmers.size()) {
	int new_kmers = num_kmers - kmers.size();

        for (int i = 0; i < new_kmers; i++) {
	    kmer = new char[kmer_len + 1];
	    kmers.push_back(kmer);
	}
    }

    for (uint i = 0; i < num_kmers; i++) {
	kmer = kmers.at(i);
        strncpy(kmer, k, kmer_len);
        kmer[kmer_len] = '\0';
        k++;
    }

    return 0;
}

int
generate_kmers(string seq, int kmer_len, int num_kmers, vector<char *> &kmers)
{
    char *kmer;
    const char *k = seq.c_str();

    for (int i = 0; i < num_kmers; i++) {
        kmer = new char[kmer_len + 1];
        strncpy(kmer, k, kmer_len);
        kmer[kmer_len] = '\0';
        kmers.push_back(kmer);
        k++;
    }

    return 0;
}


int populate_kmer_hash(string seed,int id, KmerHashMap &kmer_map, vector<char *> &kmer_map_keys, int kmer_len)
{
   //map<string, Tag *>::iterator it;
    set<string>::iterator it;
    Tag    *tag;
    vector<char *>  kmers;
    bool            exists;
  
    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
   int num_kmers = strlen(seed.c_str()) - kmer_len + 1;
   generate_kmers(seed, kmer_len, num_kmers, kmers);
        for (int j = 0; j < num_kmers; j++) {
            exists = kmer_map.count(kmers[j]) == 0 ? false : true;
         kmer_map[kmers[j]].push_back(id);   
            if (exists)
                delete [] kmers[j];
            else
                kmer_map_keys.push_back(kmers[j]);
        }   
        kmers.clear();
    return 0;
}
int 
free_kmer_hash(KmerHashMap &kmer_map, vector<char *> &kmer_map_keys) 
{
    KmerHashMap::iterator kmer_map_it;
    
    for (kmer_map_it = kmer_map.begin(); kmer_map_it != kmer_map.end(); ) 
	{
		KmerHashMap::iterator this_it = kmer_map_it++;
		this_it->second.erase(this_it->second.begin(), this_it->second.end());
  		kmer_map.erase(this_it);
  			
    }

    
    kmer_map.clear();

    for (uint i = 0; i < kmer_map_keys.size(); i++) {
        delete [] kmer_map_keys[i];
    }
    kmer_map_keys.clear();

    return 0;
}


void uclust_dist_calc(std::map<int, Tag *> &ptags, std::map<int, std::vector<int> > &tags, vector<int> &seed_list,int distance)
{

std::map<std::string, std::set<std::string> >::iterator tags_it;
std::map<std::string, std::set<std::string> > seed_dist_list;
std::set<std::string>::iterator tags_list_it;
//vector<string> seed_list;
std::map<int,int> dist_count;
vector<int>::iterator cls_tags_it;
vector<int>::iterator cls_tags_it1;
std::map<int, Tag *>::iterator ptags_it;

int d=0, seed_distance = distance*2 ;
vector<string> keys;
vector<string> t_keys;
vector<int> tags_list;
set<int>::iterator seed_matches_it;

vector<pair <string, string > > seed_pairs; 
vector<pair <string, string > >::iterator seed_pairs_it; 

string tag_1, tag_2;
Tag *tags_1, *tags_2;  

std::map<int, std::vector<int> >::iterator tags_it1;


cerr <<"\n" << "Total Number of Clusters: " << tags.size() << "\n";
cerr << "Calculationg distances between unique stacks within each cluster..." << "\n";


for ( tags_it1 = tags.begin(); tags_it1 != tags.end(); tags_it1++)
{

tags_list = tags_it1->second;
for ( int j=0; j < tags_it1->second.size(); j++)
  {
  
  tag_1 = ptags[tags_list[j]]->seq;
    for ( int k=0; k <tags_it1->second.size(); k++)
  {
    tag_2 = ptags[tags_list[k]]->seq;
     if ( j == k)
      {
        ptags[tags_list[j]]->add_dist(tags_list[k],0);
        
      }
      else
      {
      d = dist(tag_1, tag_2, distance);
      dist_count[d]++;
      if ( d != -1)
      {
     ptags[tags_list[j]]->add_dist(tags_list[k],d);
     }
      
      }
   }  


  }
   


}


d =0;



int sc=0;
cerr << "Calculationg distances between cluster seed sequences... \n";


#pragma omp parallel private(tags_1, tags_2)
 { 
  
#pragma omp for  schedule(dynamic) 
for (int l=0; l < seed_list.size()-1 ; l++)
{

 //cerr << "Calculationg distances between cluster seed sequences... " << l << "       \r";
tags_1 = ptags[seed_list[l]];

for (int m=l+1; m < seed_list.size() ; m++)
{

tags_2 = ptags[seed_list[m]];
d = dist(tags_1,tags_2, seed_distance);
if ( d > 0)
{
tags_1->seed_matches.insert(tags_2->id);
 
 

 } 
}  
  
}

}

//cout <<"\n" << "seed_pairs:  " << seed_pairs.size() << "\t" << sc << "\n";
cerr << "Calculationg distances between unique stacks among similar clusters ...\n";

int l = 0;
for ( ptags_it = ptags.begin(); ptags_it != ptags.end(); ptags_it++)
{

for ( seed_matches_it = ptags_it->second->seed_matches.begin(); seed_matches_it != ptags_it->second->seed_matches.end(); seed_matches_it++)
{
l++;
 //cerr << "Calculationg distances between unique stacks among similar clusters ... " << l << "       \r";
//seed_dist_list[ptags_it->first].insert(*seed_matches_it);
for ( cls_tags_it = tags[ptags_it->first].begin();cls_tags_it != tags[ptags_it->first].end(); cls_tags_it++)
  {
 for ( cls_tags_it1 = tags[*seed_matches_it].begin();cls_tags_it1 != tags[*seed_matches_it].end(); cls_tags_it1++)
  {
      d = dist(ptags[*cls_tags_it], ptags[*cls_tags_it1], distance);
      if ( d != -1) ptags[*cls_tags_it]->add_dist(*cls_tags_it1,d);

    

    }
  } 
 }

}



}

int calc_kmer_distance(std::map<int, Tag *> &ptags, int query, map<int, vector<int> > &clusters, KmerHashMap &kmer_map,vector<char *>  &kmer_map_keys, int kmer_len,int num_kmers,int min_hits ) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    Tag   *tag_1, *tag_2;
    int count;
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
	KmerHashMap::iterator h;
	vector<char *>        query_kmers;
	initialize_kmers(kmer_len, num_kmers, query_kmers);
    tag_1 = ptags[query];
    generate_kmers_lazily(tag_1->seq.c_str(), kmer_len, num_kmers, query_kmers);
    map<int, int> hits;
    int d;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (int j = 0; j < num_kmers; j++) {
                h = kmer_map.find(query_kmers[j]);

                if (h != kmer_map.end())
                    for (uint k = 0; k <  h->second.size(); k++)
                      {
                        hits[h->second[k]]++;
                      }
            }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {

                if (hit_it->second >= min_hits)
                {

                tag_2 = ptags[hit_it->first];

                d = dist(tag_1, tag_2, max_nwk_dist);

                if ( d != -1) 
                {
                 clusters[tag_2->id].push_back(tag_1->id);
                 break; 
                }
                
                }
                
            }
    
	for (uint j = 0; j < query_kmers.size(); j++)
	  delete [] query_kmers[j];
    return d;
}

int build_kmer_clusters (std::map<int, Tag *> &ptags) {

KmerHashMap   kmer_map;
   //vector<string> kmer_map_keys;
   vector<char *> kmer_map_keys;
    map<int, Tag *>::iterator it;
    vector<int> seeds;
    map<int, vector<int> > clusters;
     Tag *tags_1, *tags_2;  
    vector<int> keys;

    int d,i=0,l=0;;
 int con_len   = strlen(ptags.begin()->second->seq.c_str());

if ( con_len > 50)
 {
    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
   
    int kmer_len,num_kmers,min_hits;
    
    
    kmer_len = determine_kmer_length(con_len, max_nwk_dist);
    num_kmers = con_len - kmer_len + 1;
    min_hits = calc_min_kmer_matches(kmer_len, max_nwk_dist, con_len, true);
    

    
    seeds.push_back(ptags.begin()->first);
   // clusters[ptags.begin()->first].push_back(ptags.begin()->first);
    populate_kmer_hash(ptags.begin()->second->seq, ptags.begin()->first,kmer_map, kmer_map_keys, kmer_len);
    cerr << "Clustering Unique Stacks ..." << "\n";  
    for (it = ptags.begin(); it != ptags.end(); it++)  {
    i++;
    if (i % 10 == 0) 
    cerr << "Calculationg distances for stack  " << i << "       \r";
    d = calc_kmer_distance(ptags,it->first, clusters,kmer_map,kmer_map_keys,kmer_len,num_kmers,min_hits);
     if ( d == -1) 
    {
    clusters[it->first].push_back(it->first);
     populate_kmer_hash(it->second->seq,it->first, kmer_map, kmer_map_keys, kmer_len);
     seeds.push_back(it->first);
   
    }

    
    }
   free_kmer_hash(kmer_map,kmer_map_keys); 
}   
   
else
{   
 for (it = ptags.begin(); it != ptags.end(); it++)  {

seeds.push_back(it->first);
break;
}


 for (it = ptags.begin(); it != ptags.end(); it++)  {
 
 keys.push_back(it->first);
 
 }
 int j=0;
 
 for (int k = 0; k < keys.size(); k++)
 {
  j=0;
  tags_1 = ptags[keys[k]];
 l++;
if (l % 10 == 0) 
cerr << "Calculationg distances for stack  " << l << "       \r";
 

 
 
 #pragma omp parallel private(tags_2)
 { 
  
#pragma omp for  schedule(dynamic)
 for ( int i=0; i <  seeds.size(); i++)
 {
 tags_2 = ptags[seeds[i]];
 d = dist(tags_1,tags_2, max_nwk_dist);
 if ( ( d >= 0) && (j == 0) )
{

 clusters[tags_2->id].push_back(tags_1->id);
 j=1;
  //break;

 } 

 
 }
 }
 
 


if ( j == 0)
{  
clusters[tags_1->id].push_back(tags_1->id);

seeds.push_back(tags_1->id);
}

 
 }
} 
 

cout << "\n clusters: " <<  clusters.size() << "\n";   
uclust_dist_calc(ptags,clusters,seeds, max_nwk_dist);

   
    return 0;

}


std::vector <std::string> read_directory( const std::string& path = std::string() )
{
	std::vector <std::string> result;
	dirent* de;
	DIR* dp;
	//errno = 0;
	dp = opendir( path.empty() ? "." : path.c_str() );
	if (dp)
    {
		while (true)
		{
			//errno = 0;
			de = readdir( dp );
			if (de == NULL) break;
			result.push_back( std::string( de->d_name ) );
		}
		closedir( dp );
		std::sort( result.begin(), result.end() );
    }
	return result;
}
int 
write_simple_output(CTag *tag, ofstream &cat_file, ofstream &snp_file, ofstream &all_file) 
{
    vector<SNP *>::iterator           snp_it;
    map<string, int>::iterator        all_it;
    set<pair<int, int> >::iterator src_it;
    string sources;
    int batch_id=1;
    for (src_it = tag->sources.begin(); src_it != tag->sources.end(); src_it++) {
	stringstream s; 
	s << (*src_it).first << "_" << (*src_it).second << ",";
	sources += s.str();
    }
    sources = sources.substr(0, sources.length() - 1);
    if ( sources == "") return 0;
    cat_file << 
	"0"          << "\t" << 
	batch_id     << "\t" <<
	tag->id     << "\t" <<
     "\t" <<
    "0" << "\t" <<
    "+" << "\t" <<
	"consensus"  << "\t" <<
	"0"          << "\t" <<
	sources      << "\t" <<
	tag->con     << "\t" << 
        0            << "\t" << 
        0            << "\t" <<
        0            << "\t" <<
        0            << "\n";

    //
    // Output the SNPs associated with the catalog tag
    //
    for (snp_it = tag->snps.begin(); snp_it != tag->snps.end(); snp_it++) {
   if ((*snp_it)->type == snp_type_het)
    {
	snp_file << "0"    << "\t" << 
	    batch_id       << "\t" <<
	    tag->id        << "\t" << 
	    (*snp_it)->col << "\t";
    snp_file << "E\t";    
	
	snp_file << 
	    (*snp_it)->lratio << "\t" << 
	    (*snp_it)->rank_1 << "\t" << 
	    (*snp_it)->rank_2 << "\t" << 
	    ((*snp_it)->rank_3 == 0 ? '-' : (*snp_it)->rank_3) << "\t" << 
	    ((*snp_it)->rank_4 == 0 ? '-' : (*snp_it)->rank_4) << "\n";
	 }  
	 delete *snp_it;
    }
   
    //
    // Output the alleles associated with the two matched tags
    //
    for (all_it = tag->alleles.begin(); all_it != tag->alleles.end(); all_it++)
    {
	all_file << 
	    "0"           << "\t" << 
	    batch_id      << "\t" <<
	    tag->id      << "\t" <<
	    all_it->first.c_str() << "\t" <<
            "0"           << "\t" <<    
            "0"           << "\n";     
   }
    return 0;
}

int write_catalog(map<int, CTag *> &catalog, string path) {
    map<int, CTag *>::iterator i;
    CTag  *tag;
    set<int> matches;

   // bool gzip = (in_file_type == FileT::gzsql) ? true : false;
 bool gzip = false;
    //
    // Parse the input file names to create the output file
    //
    stringstream prefix; 
    prefix << path << "batch_1";

    string tag_file = prefix.str() + ".catalog.tags.tsv";
    string snp_file = prefix.str() + ".catalog.snps.tsv";
    string all_file = prefix.str() + ".catalog.alleles.tsv";

    if (gzip) {
        tag_file += ".gz";
	snp_file += ".gz";
	all_file += ".gz";
    }

    //
    // Open the output files for writing.
    //
    gzFile   gz_tags, gz_snps, gz_alle;
    ofstream tags, snps, alle;
    if (gzip) {
        gz_tags = gzopen(tag_file.c_str(), "wb");
	if (!gz_tags) {
	    cerr << "Error: Unable to open gzipped catalog tag file '" << tag_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
        #if ZLIB_VERNUM >= 0x1240
	gzbuffer(gz_tags, libz_buffer_size);
	#endif
	gz_snps = gzopen(snp_file.c_str(), "wb");
	if (!gz_snps) {
	    cerr << "Error: Unable to open gzipped catalog snps file '" << snp_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
        #if ZLIB_VERNUM >= 0x1240
	gzbuffer(gz_snps, libz_buffer_size);
	#endif
	gz_alle = gzopen(all_file.c_str(), "wb");
	if (!gz_alle) {
	    cerr << "Error: Unable to open gzipped catalog alleles file '" << all_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
        #if ZLIB_VERNUM >= 0x1240
	gzbuffer(gz_alle, libz_buffer_size);
	#endif
    } else {
	tags.open(tag_file.c_str());
	if (tags.fail()) {
	    cerr << "Error: Unable to open catalog tag file for writing.\n";
	    exit(1);
	}
	snps.open(snp_file.c_str());
	if (snps.fail()) {
	    cerr << "Error: Unable to open catalog SNPs file for writing.\n";
	    exit(1);
	}
	alle.open(all_file.c_str());
	if (alle.fail()) {
	    cerr << "Error: Unable to open catalog alleles file for writing.\n";
	    exit(1);
	}
    }

    //
    // Record the version of Stacks used and the date generated as a comment in the catalog.
    //
    // Obtain the current date.
    //
    stringstream log;
    time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
   log << "# RAD2POP version 1.0; catalog generated on " << date << "\n"; 
    if (gzip) {
        gzputs(gz_tags, log.str().c_str());
        gzputs(gz_snps, log.str().c_str());
        gzputs(gz_alle, log.str().c_str());
    } else {
        tags << log.str();
	snps << log.str();
	alle << log.str();
    }

    for (i = catalog.begin(); i != catalog.end(); i++) {
	tag = i->second;

	//if (gzip)
	   // write_gzip_output(tag, gz_tags, gz_snps, gz_alle);
	//else 
	 if (tag->con != "")
	    write_simple_output(tag, tags, snps, alle);
	    
	 delete i->second;  
    }

    if (gzip) {
	gzclose(gz_tags);
	gzclose(gz_snps);
	gzclose(gz_alle);
    } else {
	tags.close();
	snps.close();
	alle.close();
    }

    return 0;
}

int mergeStacks(std::map<int, Tag *>  &ptags,vector<string> &samples, int dist, int cov,  ofstream &stats_fh)
{

std::map<int, CTag *> ctags;
std::map<int,set<int> > ctag_alleles_map;
std::vector<int>::iterator vec_it;
 std::map<int, std::vector<int> >::iterator it;  
std::stringstream ss;
ss << dist;
string dist_str = ss.str();

std::stringstream ss1;
ss1 << cov;
string cov_str = ss1.str();

string path = out_path+"M"+dist_str+"_m"+cov_str;
string param = "M"+dist_str+"_m"+cov_str;
string cmd= "mkdir -p " + path;
int ret=system(cmd.c_str());
path = path+"/" ;

string ind_id = "";
int pos, pos1;	 

std::map<int, std::vector<int> >::iterator merged_it;	 
std::map<int, Tag *>::iterator ptags_it;
vector<pair<int, int> >::iterator it2;
map<int, set<int> >::iterator  matches_it;
set<int>::iterator matches_val_it;
map<string, int> pop_samples;

int size =0, merged_sec_cnt=0;
int li =-1;

for(int sample=1; sample <= samples.size(); sample++)
{
pop_samples[samples[sample-1].substr(0,3)]++;

std::map<int, std::vector<int> > merged;
map<int, set<int> >  matches;
std::vector<int> merged_vec;

for (ptags_it = ptags.begin(); ptags_it != ptags.end(); ptags_it++) {

                pos=0;
                for( vec_it = ptags_it->second->sample_ids.begin(); vec_it != ptags_it->second->sample_ids.end(); vec_it++)
                {
                   if (sample== *vec_it) break;
                   pos++;
                }
                
if ( ptags_it->second->sample_ids.size() != pos)
{
for (it2 = ptags_it->second->dist.begin(); it2 != ptags_it->second->dist.end(); it2++)
{  
  pos1 = 0;
  for( vec_it = ptags[it2->first]->sample_ids.begin(); vec_it != ptags[it2->first]->sample_ids.end(); vec_it++)
                {
                   if (sample== *vec_it) break;
                   pos1++;
                }
  if ( ptags[it2->first]->sample_ids.size() != pos1) 
  {
    if  ( ( ptags_it->second->cov.at(pos) >= cov) && ( it2->second <= dist) )
				{   
				       if (ptags[it2->first]->cov.at(pos1) >= cov)
				      {
					merged[ptags_it->first].push_back(it2->first);
					if (ptags_it->first != it2->first) 
					{
				    merged[it2->first].push_back(ptags_it->first);
				    }	
				   }	
					
				}
   				
  
  else if  ( ( ptags_it->second->cov.at(pos) < cov) && ( it2->second <= dist+2)) 	
                {
                  if ((ptags_it->second->dist.size() > 1 ) && (ptags[it2->first]->cov.at(pos1) >= cov))
                  {
                    matches[ptags_it->first].insert(it2->first);
                    
                   } 
                    
  			    }
 } 			    
  			    
 }
 }
 
 }
size =0;
li = -1;
merged_sec_cnt=0;
vector<int> mkeys;

		for (merged_it = merged.begin(); merged_it != merged.end(); merged_it++) 
		{
		
			mkeys.push_back(merged_it->first);
		}
		
		for (int k=0; k < mkeys.size(); k++)
		{   
			
			size =  merged[mkeys[k]].size();
			for ( int i =0; i < size; i++)
			{
				
				li = merged[mkeys[k]].at(i);
				if ( mkeys[k] == li) continue;
				if (merged[li].size() >= 1) 
				{   
					merged[mkeys[k]].insert( merged[mkeys[k]].end(), merged[li].begin(), merged[li].end());
					
					merged[li].clear();
					std::sort( merged[mkeys[k]].begin(), merged[mkeys[k]].end() );
					merged[mkeys[k]].erase( std::unique( merged[mkeys[k]].begin(), merged[mkeys[k]].end() ), merged[mkeys[k]].end() );
					size =  merged[mkeys[k]].size();
					i = 0;
				}
			} 
		
		}
	
 

     
   for (merged_it = merged.begin(); merged_it != merged.end(); merged_it++) 
       {
       if (merged_it ->second.size() > max_stacks)
        merged[merged_it->first].clear();
       
        
       }
 
  for (matches_it=matches.begin(); matches_it!=matches.end(); matches_it++)
       {
       
       for ( matches_val_it = matches_it->second.begin(); matches_val_it != matches_it->second.end(); matches_val_it++)
         {
          if ( merged[*matches_val_it].size() > 0)
          {
             merged[*matches_val_it].push_back(matches_it->first);
             merged_sec_cnt++;
             break;
          }
         
         }
       
       
       }

       cerr  << "# of sec merged stacks" << "\t"<< merged_sec_cnt << '\n';
		ind_id = samples[sample-1];
		call_consensus(merged, ptags,ctags,ctag_alleles_map,sample, ind_id, path, cov, param, stats_fh);
        cerr << "catalog size: " << ctags.size() << "\n";
        merged.clear();
        matches.clear();
}

write_catalog(ctags, path);
ctags.clear();
return 0;
}



int main (int argc, char* argv[]) 
{
 parse_command_line(argc, argv);
 

  std::cerr  << " Parameter sweep mode: " << (psweep == true ? "on"  : "off") << "\n"
	          << " Maximum distance (in nucleotides) allowed between stacks to form network: "<< max_nwk_dist  << "\n"
	          << " Minimum coverage depth: " << min_merge_cov << "\n" 
              << " Maximum number stacks per locus: " << max_stacks << "\n"
              << " Minimum sample percentage: " << min_sam_per << "\n"
              << " Minimum average coverage depth: " << min_depth << "\n";
 
 
 #ifdef _OPENMP
 omp_set_num_threads(num_threads);
 #endif
  
  
 map<int,int> sample_count;
 map<int,int>::iterator sample_count_it;
   

 map<string,int>::iterator radtags_it;
 std::map<string, Tag *> ptags;
 std::map<int, Tag *> tags;
 std::map<string, Tag *>::iterator ptags_it;
  std::map<int, Tag *>::iterator tags_it;
 Tag *tag;
 DNASeq *seq;
 std::vector <std::string>  files;
 std::vector <std::string>::iterator files_it;
 int i =0, sample_id = 1;
 vector<string> samples;
 vector<int>::iterator cov_it;
 string stats_file; 
 int old_id =0; 
 files = read_directory (in_file);
 
  time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo); 
    stats_file = out_path + "De_novo_stats";
    ofstream stats_fh(stats_file.c_str(), ofstream::out);
     
     if (stats_fh.fail()) {
        cerr << "Error opening Stats file '" << stats_file << "'\n";
	exit(1);
    }
	
	stats_fh << "# RAD2POP version 1.0 generated on " << date << "\n";
    stats_fh << "Sample ID" << "\t" << "Parameters" << "\t" << "Number of Loci " << "\t" << "Number of Polymorphic Loci " << "\t" << "Number of SNPs " << "\n";

 
 
for ( files_it = files.begin(); files_it != files.end(); files_it++)
	{   
  		if ( (files_it->find("fq") != std::string::npos) || (files_it->find("fa") != std::string::npos) || (files_it->find("gz") != std::string::npos))
			{  
    			i=0;
    			//old_id =0; 
    			
    			map<string,int> radtags;
    			string filename = in_file+"/"+ *files_it;
    			string samplename = *files_it;
    			size_t lastindex = samplename.find_last_of("."); 
    			string rawname = samplename.substr(0, lastindex);
    			string popname = rawname.substr(0,3); 
    			cerr <<  "Identifying unique stacks for:  " << samplename << "\n";
    			samples.push_back(rawname);   
				load_radtags(radtags,filename); 
    
	
 			for (radtags_it = radtags.begin(); radtags_it != radtags.end(); radtags_it++) 
 				{
       				if (i % 10 == 0) cerr << "Loading stack " << i << "       \r";
    				if (radtags_it->second > 1) 
    					{
             				ptags_it = ptags.find(radtags_it->first);
             				if (ptags_it == ptags.end())
             					{ 
             						tag = new Tag;
             						tag -> id = old_id;
             						tag -> seq= radtags_it->first;
             						tag -> len= tag -> seq.length(); 
             						tag -> cov.push_back(radtags_it->second);
             						tag -> sample_ids.push_back(sample_id);
             						tag -> pop_ids.insert(popname);;
             						ptags[radtags_it->first] = tag; 
             						old_id = tag -> id+1;  
             					}
             				else
             					{
             						tag = ptags[radtags_it->first];
             						tag -> cov.push_back(radtags_it->second);  
             						tag -> pop_ids.insert(popname); 
             						tag -> sample_ids.push_back(sample_id);
             					}
	         				i++;     
           				}      
				} 
			}
			else continue;
	sample_id++;
}



int avg_cov = 0, sum_cov = 0,sample_size = 0, min_sam_fil = (sample_id-1)*min_sam_per;


cerr << "Number of Unique stacks (before filter):" << ptags.size() <<"\n"; 
cerr << "Minimum Samples: " <<  min_sam_fil <<"\n"; 
cerr << "Minimum Average Coverage: " <<  min_depth << "\n"; 

if( ptags.size() > 0)
{
for (ptags_it = ptags.begin(); ptags_it != ptags.end(); ) 
	{
		std::map<string, Tag *>::iterator this_it = ptags_it++;
		sample_size = this_it->second->sample_ids.size();
		for ( cov_it = this_it->second->cov.begin(); cov_it != this_it->second->cov.end(); cov_it++)
			{
				sum_cov = sum_cov + *cov_it;
			}
		avg_cov = sum_cov/sample_size;
		sum_cov = 0;
	
	    	
		if (this_it->second->sample_ids.size() <= min_sam_fil ) 
			{
 				if (avg_cov <= min_depth) 
  				 {      
  				        delete this_it->second; 
  						ptags.erase(this_it);
  				 }
			}

	}

for (ptags_it = ptags.begin(); ptags_it != ptags.end(); ptags_it++) 
	{
	 tags[ptags_it->second->id] =ptags_it->second;
	
	
	}

ptags.clear();


cerr << "Number of Unique stacks (after filter):" << tags.size() <<"\n";

build_kmer_clusters(tags);
 
if ( psweep)
{ 
for ( int i = 2; i <= max_nwk_dist-2; i++)
	{
		for ( int j=3; j <= min_merge_cov; j++)
			{
				cerr << "\n"<< "Maximum Nucleotide distance M:" << "\t" <<i;
				cerr << "\n"<< "Minimum Coverage m:" << "\t" << j;
				mergeStacks(tags, samples,i,j,stats_fh);
			}
	}
}
else
{
mergeStacks(tags, samples,max_nwk_dist-2,min_merge_cov, stats_fh);
}
cerr << "\n"<< "Process Completed !!!" << "\n";
}

else 
{
cerr << "\n"<< "Process Aborted!! No RADTags Loaded !!!" << "\n";
return 0;
}
return 0;

}


int load_radtags(map<string,int> &radtags, string filename) {
    
    Input *fh;
   if (in_file_type == "fasta")
        fh = new Fasta(filename.c_str());
    else if (in_file_type == "fastq")
        fh = new Fastq(filename.c_str());
    else if (in_file_type == "gzfasta")
        fh = new GzFasta(filename.c_str());
    else if (in_file_type == "gzfastq")
        fh = new GzFastq(filename.c_str());

   // cerr << "Parsing " << filename.c_str() << "\n";
    long  int corrected = 0;
    long  int i         = 0;
    short int seql      = 0;
    short int prev_seql = 0;
    bool  len_mismatch  = false;

    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {

	prev_seql = seql;
	seql      = 0;
    string seq="";
	for (char *p = c.seq; *p != '\0'; p++, seql++)
	    switch (*p) {
	    case 'N':
	    case 'n':
	    case '.':
		*p = 'A';
		corrected++;
		
	    }
	if (seql != prev_seql && prev_seql > 0) len_mismatch = true; 
    radtags[string(c.seq)]++;
      i++;
    }
 

    if (len_mismatch)
	cerr << "  Warning: different sequence lengths detected.\n";

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}



