//
//  main.cpp
//  sthlm
//
//  Created by Jilin Zhang on 2/14/17.
//  Copyright © 2017 Jilin Zhang. All rights reserved.
//



#include <iostream>
#include <fstream>
#include <string>
#include "unordered_map"
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <boost/program_options.hpp>
//#include "omp.h"
#include "kseq.h"
#include <zlib.h>

KSEQ_INIT(gzFile,gzread);

//
struct Unitscore{
    long count;
    double discore;
    double monoscore;
};

std::unordered_map<std::string, Unitscore > kstore;


//g++ -std=c++11 main.cpp -o rmap
//g++ -std=c++11 main.cpp -o rmap -lboost_program_options -lz
/*
 basemap[] works by storing a very small array that maps a base to
 its complement, by dereferencing the array with the ASCII char's
 decimal value as the index
 (int) 'A' = 65;
 (int) 'C' = 67;
 (int) 'G' = 71;
 (int) 'T' = 84;
 (int) 'a' = 97;
 (int) 'c' = 99;
 (int) 'g' = 103;
 (int) 't' = 116;
 (int) 'N' = 78;
 (int) 'U' = 85;
 (int) 'u' = 117;
 for example: basemap['A'] => basemap[65] => 'T' etc.
 */
static const char basemap[255] =
{
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
    '\0', '\0', '\0', '\0', '\0',  'U', '\0',  'G', '\0', '\0', /*  60 -  69 */
    '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
    '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0',  'u', '\0',  'g', /*  90 -  99 */
    '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
    '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
    '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
};

void fetch_matrix(std::unordered_map<std::string, std::vector<double> > *map, std::string str, int ind, double count);
void show_matrix(std::unordered_map<std::string, std::vector<double> > *x);
void kmerizing (std::string *id, std::string* seq,  int kmer, int gap, std::string* name, std::unordered_map<std::string, std::vector<double> > *map);
double monoscoring(std::string subseq, std::unordered_map<std::string, std::vector<double> > *map, int gap_len=0, int gap_pos=0);
double discoring(std::string subseq, std::unordered_map<std::string, std::vector<double> > *map, int gap_len=0, int gap_pos=0);
inline double singleStringScore(double x, std::string* s, long p, std::unordered_map<std::string, std::vector<double> > *m);
inline double pairStringScore(double x, std::string* s, long sp, std::unordered_map<std::string,std::vector<double> > *m);
void Matrixfile(const char*x, std::unordered_map<std::string, std::vector<double> > *y, double pseudo_count);
int listParser(const char* list_name, std::string* names, int* length, std::unordered_map<std::string, std::vector<double> >* map,double pcount);
std::string seq_revcomp(std::string str);
void GenomeMap(std::string *id, std::string *seq, long monoLen, std::string *mname, std::unordered_map<std::string, std::vector<double> >* map );
void hashOut();
void unCenterCorr();



template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
         copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
         return os;
}

int main(int argc,   char * argv[]) {
    
    std::unordered_map<std::string, std::vector<double> > matrixmap[100];
    int motif_len[100];
    std::string names_motif[100];
    double pseudo_count=0.00001;
    int motif_counts=0;
    char list_name[100],fasta_name[100];
    int kmer_length, gap_length;
    bool kmerflag=false ,outflag=false,corrflag=false, revcomflag=false;
    
    
    boost::program_options::options_description desc("Allowed options:");
    desc.add_options()
    ("help", "produce help message")
    ("fasta", boost::program_options::value<std::string>(), "input sequence file")
    ("list", boost::program_options::value<std::string> (), "specify a list of motifs")
    ("motif", boost::program_options::value<std::string>(), "specify one motif file")
    ("kmer,k", boost::program_options::value<int>(&kmer_length)->default_value(8), "specify the kmer length")
    ("gap,g", boost::program_options::value<int>(&gap_length)->default_value(0), "specify the gap length")
    ("kmerflag,f", boost::program_options::bool_switch(&kmerflag), "flag to switch to the kmer counting-score procedure")
    ("revcomp,r",boost::program_options::bool_switch(&revcomflag),"reverse compliment sequence")
    ("output,o",boost::program_options::bool_switch(&outflag),"output the kmer count-score")
    ("corr,c",boost::program_options::bool_switch(&corrflag),"output the correlation score (un centered cosine)")
    ;
    boost::program_options::variables_map vm;
    boost::program_options::positional_options_description p;
    p.add("fasta",-1);
    boost::program_options::store(boost::program_options::command_line_parser(argc,argv).options(desc).positional(p).run(), vm);
    boost::program_options::notify(vm);
    
    if (vm.count("list"))
    {
        const char* temx=vm["list"].as<std::string >().c_str();
        motif_counts= listParser(temx,names_motif, motif_len, matrixmap, pseudo_count);
    }
    if (vm.count("motif")) {
        const char* temname =vm["motif"].as<std::string>().c_str();
        motif_counts=1;
        strcpy(list_name,temname);
        names_motif[0]=(std::string) list_name;
        Matrixfile(temname, matrixmap, pseudo_count);
    }
    if (vm.count("fasta"))
    {
        const char * fastaname=vm["fasta"].as<std::string>().c_str();
        motif_len[0]=matrixmap[0].at("A").size();
        strcpy(fasta_name, fastaname );
    }
    if (vm.count("help")) {
        std::cout << "Usage:" << argv[0] << " [options]\n";
        std::cout << desc;
        return 0;
    }

    

// Load fasta sequence one by one;
    gzFile fp;
    kseq_t *seq;
    int l;
    fp=gzopen(fasta_name, "r");
    seq=kseq_init(fp);
    while((l=kseq_read(seq)) >=0){
        std::string seqnametmp, seqtmp;
        seqnametmp= &seq->name.s[0];
        seqtmp= &seq->seq.s[0];
        if(revcomflag)
            seqtmp=seq_revcomp(seqtmp);
        
        for(int i=0;i < motif_counts ; i++){
            if(kmerflag)
                kmerizing(&seqnametmp, &seqtmp,kmer_length,gap_length,&names_motif[i], &matrixmap[i]);
            else
                GenomeMap(&seqnametmp, &seqtmp,motif_len[i],&names_motif[i], &matrixmap[i]);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    if(kmerflag && outflag)
        hashOut();
    if(kmerflag && corrflag)
        unCenterCorr();

    return 0;
}

//output the kmer count-score table
void hashOut(){
    for(auto it=kstore.begin(); it !=kstore.end(); it++)
        std::cout<<it->first << "\t" << it->second.count<< "\t" << it->second.discore <<std::endl;
//        std::cout<<it->first << "\t" << it->second.count<< "\t" << it->second.discore << "\t" << it->second.monoscore <<std::endl;

}

// function to calculate the uncentred the correlation, cosine
void unCenterCorr(){
    double sumx,sumy, sumz, sumxz, sumxy;
    for(auto it=kstore.begin(); it !=kstore.end(); it++){
        sumx+=it->second.count * it->second.count;
        sumy+=it->second.discore * it->second.discore ;
        sumz+=it->second.monoscore * it->second.monoscore  ;
        sumxy+= it->second.count * it->second.discore;
        sumxz+= it->second.count * it->second.monoscore;
    }
//    std::cout << sumxy/sqrt(sumx *sumy)<< "\t" << sumxz/sqrt(sumx *sumz) << std::endl;
    std::cout << sumxy/sqrt(sumx *sumy)<< std::endl;
}

// function to score all the locus in the sequence
void GenomeMap(std::string *id, std::string *seq, long monoLen, std::string *mname, std::unordered_map<std::string, std::vector<double> >* map ){
    unsigned long slen =seq->size();
    
    std::transform(seq->begin(), seq->end(), seq->begin(),::toupper);
    std::replace(seq->begin(),seq->end(),'T','U');
    
    for(long i=0; i < slen - monoLen +1 ; i++){
        std::string subseq=seq->substr(i,monoLen);
        if(subseq.find("N") != std::string::npos)
            continue;
        double mscore=0,dscore=0,rc_mscore=0,rc_dscore=0;
        std::string rc_subseq=seq_revcomp(subseq);
        
        mscore=singleStringScore(mscore, &subseq, 0, map);
        dscore= pairStringScore(dscore, &subseq,0, map);
        rc_mscore=singleStringScore(rc_mscore, &rc_subseq, 0,map);
        rc_dscore=pairStringScore(rc_dscore,&rc_subseq,0,map);
        std::cout << *id << "\t"<< i << "\t" << std::setprecision(5) <<mscore << std::fixed <<":"  << std::setprecision(5) <<dscore <<std::fixed << "\t+\t" << *mname <<std::endl;
        std::cout << *id << "\t"<< i << "\t" << std::setprecision(5) <<rc_mscore <<std::fixed << ":"  << std::setprecision(5) <<rc_dscore << std::fixed <<"\t-\t" <<*mname <<std::endl;
    }    
}

// read the motif list
int listParser(const char* list_name, std::string* names, int* length, std::unordered_map<std::string, std::vector<double> >* map,double pcount){
    std::fstream motif_list;
    char motif_file_name[100];
    motif_list.open(list_name);
    int motif_line=0;
    
    while( motif_list.getline( motif_file_name,100) ){
        Matrixfile(motif_file_name, &map[motif_line],  pcount);
//        std::cout << motif_file_name << std::endl;
        length[motif_line]=(int) (map[motif_line])["A"].size();
        names[motif_line]=(std::string) motif_file_name;
        ++motif_line;
    }
//    std::cout << motif_line << " motifs loaded \n"<<std::endl;
    return motif_line;
}

//process the matrix input
void Matrixfile(const char* x, std::unordered_map<std::string, std::vector<double> > *y, double pseudo_count){
    std::fstream file;
    file.open(x);
    
    while(!file.is_open()){
        std::cout << "Something wrong with the filename" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    int line_count=0;
    std::string iss;
    while(getline(file, iss)){
        line_count++;
        if(line_count <=2){
            continue;
        }
        else if(line_count >=3  and line_count <=18){
            fetch_matrix(y, iss, 1, pseudo_count);
        }
        else if(line_count >=19 && line_count<=22){
            fetch_matrix(y, iss, 2, pseudo_count);
        }
    }
}

//show the matrix and output or the screen
void show_matrix(std::unordered_map<std::string, std::vector<double> > *x){
    for(auto it=(*x).begin(); it !=(*x).end(); it++){
        std::cout << it->first << "\t";
        std::vector<double> tem= it->second;
        for( int i=0; i< tem.size(); ++i)
            std::cout <<  "\t" <<  tem[i];
        std::cout << std::endl;
    }
}
std::string seq_revcomp(std::string str){
    std::string rc_seq;
    for(int i=str.size()-1; i>=0; --i){
        rc_seq += (char) basemap[(int)str[i]];
    }
    return rc_seq;
}

//store the matrix into map table
void fetch_matrix(std::unordered_map<std::string, std::vector<double> > *map, std::string str, int ind, double count){
    std::vector<std::string> matrix;
    std::stringstream ms(str);
    std::string units;
    while(ms >> units){
        matrix.push_back(units);
    }
    unsigned long matrix_length_raw=matrix.size();
    std::string base_name= matrix.back();
    long matrix_length_new;
    if(ind ==1){
        base_name=base_name.substr(0,2);
        matrix_length_new=matrix_length_raw-2;
    }
    else if(ind==2){
        base_name=base_name.substr(14,1);
        matrix_length_new=matrix_length_raw-1;
    }
    double new_matrix[matrix_length_new];
    std::vector<double> hallo (new_matrix,new_matrix+ matrix_length_new);
    for(int i=0; i< matrix_length_new; i++){
        std::string::size_type sz;
        hallo[i]= std::stod(matrix[i], &sz);
        if(ind ==1)
            hallo[i] = std::log2((hallo[i]+ count)/(0.0625 + count));
        else
            hallo[i] = std::log2((hallo[i]+ count)/(0.25+ count));
    }
    map->insert({base_name,hallo});
//    return map;
}

//function to generate the kmer count-score table
void kmerizing (std::string* id, std::string* seq,  int kmer, int gap, std::string* name, std::unordered_map<std::string, std::vector<double> > *map){
    std::transform(seq->begin(), seq->end(), seq->begin(),::toupper);
    std::replace(seq->begin(),seq->end(),'T','U');
    unsigned long seq_length=seq->size();
    int raw_kseq_len=kmer+gap;

    for(unsigned long i=0; i < seq_length-raw_kseq_len +1; i++){
        std::string subseq=seq->substr(i,raw_kseq_len);
        if(subseq.find("N") != std::string::npos)
            continue;
        
        if(gap >0){
            for (int ii=1; ii <kmer ;ii++){
                double mscore,dscore;
                mscore=monoscoring(subseq, map, gap,ii);
                dscore=discoring(subseq, map, gap,ii);
                std::string tem=subseq;
                tem.replace(ii,gap, gap, '-' );

                if(kstore.find(tem) != kstore.end()){
                    kstore[tem].count +=1;
                }
                else
                    kstore.insert({tem,{1,dscore,mscore}});
            }
        }
        else{
            double mscore,dscore,rc_mscore, rc_dscore;
            //         std::string rc_subseq=seq_revcomp(subseq);
            mscore=monoscoring(subseq, map, 0,0);
            dscore=discoring(subseq, map, 0,0);
            std::string tem=subseq;
            
            if(kstore.find(tem) != kstore.end()){
                kstore[tem].count +=1;
            }
            else
                kstore.insert({tem,{1,dscore,mscore}});
        }


//        double mscore,dscore,rc_mscore,rc_dscore;
//        std::string rc_subseq=seq_revcomp(subseq);
//        mscore=monoscoring(subseq, map, 0,0);
//        dscore=discoring(subseq, map, 0,0);
//
//        rc_mscore=monoscoring(rc_subseq, map, 0,0);
//        rc_dscore=discoring(rc_subseq, map, 0,0);
//        std::cout << *id << "\t"<< i << "\t" << std::setprecision(5) <<mscore << std::fixed <<":"  << std::setprecision(5) <<dscore <<std::fixed << "\t+\t" << *name <<std::endl;
//        std::cout << *id << "\t"<< i << "\t" << std::setprecision(5) <<rc_mscore <<std::fixed << ":"  << std::setprecision(5) <<rc_dscore << std::fixed <<"\t-\t" <<*name <<std::endl;
//
    }
}

// The scoring function using the linear model only, which has considered all the conditions
double monoscoring(std::string subseq, std::unordered_map<std::string, std::vector<double> > *map, int gap_len, int gap_pos){
    long monoLen=(*map)["A"].size();
    if(gap_len >0)
        subseq.replace(gap_pos,gap_len, gap_len, '-' );
    
    long kseq_prime_len= subseq.size();
    long cycle=kseq_prime_len+ monoLen;
    std::string valid_seq;
    double maximum_score=-1000;
    for(int p=0; p < cycle -1 ; p++){
        long pos;
        double score=0;
        if(p <kseq_prime_len-1){
            if(monoLen >= kseq_prime_len){
                valid_seq=subseq.substr(kseq_prime_len-(p+1),p+1);
                pos=0;
            }
            else if(monoLen < kseq_prime_len){
                pos=0;
                if(p <monoLen -1)
                    valid_seq=subseq.substr(kseq_prime_len-(p+1));
                else
                    valid_seq=subseq.substr(kseq_prime_len-(p+1),monoLen);
            }
        }
        else if(p >=kseq_prime_len -1){
            if(monoLen >=kseq_prime_len){
                pos=p+1-kseq_prime_len;
                if(p <monoLen-1)
                    valid_seq=subseq;
                else
                    valid_seq=subseq.substr(0,kseq_prime_len+monoLen-p-1);
                }
            else{
                valid_seq=subseq.substr(0,kseq_prime_len-(p-monoLen) -1);
                pos=p+1-kseq_prime_len;
            }
        }
        score= singleStringScore(score, &valid_seq, pos, map);
        maximum_score=(score>maximum_score)?score:maximum_score;
    }
    return maximum_score ;
}
// The scoring function for the dependency matrix
double discoring(std::string subseq, std::unordered_map<std::string, std::vector<double> > *map, int gap_len, int gap_pos){
    if(gap_len >0)
        subseq.replace(gap_pos,gap_len,gap_len,'-');
    long monoLen=(*map)["A"].size();
    long kseq_prime_len=subseq.size();
    double maximum_score=-1000;
    if(monoLen %2 ==0){ //matrix length is even
        for(int p=0; p < monoLen + kseq_prime_len -1; p++){
            std::string subs;
            double score=0;
            long pos;
            if( p< monoLen/2){
                if( p < kseq_prime_len-1){
                    subs=subseq.substr(kseq_prime_len-1-p,p+1); //for kmers longer than the 1/2 matrix length
                    pos=0;
                }
                else if(p >=kseq_prime_len-1){  //for kmers shorter than the 1/2 marix length
                    subs=subseq;
                    pos=p-kseq_prime_len+1;
                }
                score= singleStringScore(score, &subs, pos, map);
                maximum_score=(score > maximum_score)?score:maximum_score;
            }
            else if(p >= monoLen/2   && (p - monoLen/2 < kseq_prime_len -1)){
                std::string left_seq; // leftover from the kmer which has been cut for dinucleotide matrix
                long single_pos;
                if(kseq_prime_len <= monoLen && kseq_prime_len >= 0.5 * monoLen){
                    if(2*(p-(monoLen/2 -1)) <=kseq_prime_len){
                        subs=subseq.substr(kseq_prime_len-2*(p-(monoLen/2 -1)),2*(p-(monoLen/2 -1)));
                        single_pos=p+1-subs.size();
                        if(p+1 >= kseq_prime_len){
                            pos=p+1-kseq_prime_len;
                            left_seq=subseq.substr(0,kseq_prime_len-2*(p-(monoLen/2-1))); //if($kseq_prime_len-2*($p-($mono_len/2 -1)) > 0);
                        }
                        else{
                            pos=0;
                            left_seq=subseq.substr(kseq_prime_len-(p+1),p+1-2*(p-(monoLen/2-1)));
                        }
                    }
                    else{
                        subs=subseq.substr(0,2*(kseq_prime_len-(p+1-(monoLen-1)/2 )+1) );
                        pos=p-kseq_prime_len+1+subs.size();
                        single_pos=p+1-kseq_prime_len;
                        if(p >= monoLen  ){
                            std::string left_seq_pre=subseq.substr(2*(kseq_prime_len-(p-(monoLen/2-1))));
                            left_seq=left_seq_pre.substr(0,kseq_prime_len-(p+1-monoLen)-subs.size());
                        }
                        else
                            left_seq=subseq.substr(2*(kseq_prime_len-(p-(monoLen/2 -1))));
                    }
//                    std::cout << subseq<< "\t" << subs << "\t" << single_pos << "\t" << pos << std::endl;
                }
                else if(kseq_prime_len > monoLen){
                    std::string left_seq_pre;
                    if(p < monoLen-1){
                        subs=subseq.substr(kseq_prime_len-2*(p-(monoLen/2 -1)),2*(p-(monoLen/2 -1)));
                        if(2*(p-(monoLen/2 -1)) <kseq_prime_len)
                            left_seq_pre =subseq.substr(0,kseq_prime_len-2*(p-(monoLen/2 -1)));
                        left_seq=left_seq_pre.substr(kseq_prime_len - (p+1));
                        pos=0;
                        single_pos=p+1-subs.size();
                    }
                    else if( p >= monoLen-1 && p <= kseq_prime_len-1){
                        subs=subseq.substr(kseq_prime_len-(p+1),monoLen);
                        pos=0;
                        single_pos=0;
                        
                    }
                    else if(p > kseq_prime_len -1){
                        subs=subseq.substr(0,2*(kseq_prime_len-(p-(monoLen/2 -1))));
                        left_seq_pre=subseq.substr(2*(kseq_prime_len-(p-(monoLen/2 -1))));
                        left_seq=left_seq_pre.substr(0, left_seq_pre.size()-(p+1-monoLen));
                        pos=p-kseq_prime_len+1+subs.size();
                        single_pos=p+1-kseq_prime_len;
                    }
//                    std::cout << subseq<< "\t" << subs << "\t" << single_pos << "\t" << pos << std::endl;
                }
                else if(kseq_prime_len <  0.5 * monoLen){
                    if(2*(p-(monoLen/2 -1)) <=kseq_prime_len){
                        subs=subseq.substr(kseq_prime_len-2*(p-(monoLen/2 -1)),2*(p-(monoLen/2 -1)));
                        pos=p+1-kseq_prime_len;
                        left_seq=subseq.substr(0,kseq_prime_len-2*(p-(monoLen/2 -1)));
                        single_pos=p+1-subs.size();

                    }
                    else{
                        subs=subseq.substr(0,2*(kseq_prime_len-(p-(monoLen/2 -1))));
                        pos=p-kseq_prime_len+1+subs.size();
                        left_seq=subseq.substr(2*(kseq_prime_len-(p-(monoLen/2 -1))));
                        single_pos=p+1-kseq_prime_len;
                    }
                }

                score=pairStringScore(score,&subs,single_pos,map);
//                std::cout << subseq<< "\t" << subs << "\t" << single_pos << "\t" << pos << std::endl;

                if(!left_seq.empty())
                    score= singleStringScore(score, &left_seq, pos, map);
                maximum_score=(score > maximum_score)?score:maximum_score;
            }
            else if(p >= monoLen/2 + kseq_prime_len -1){
                if(p <monoLen){
                    subs=subseq;
                    pos=p+1-kseq_prime_len;
                }
                else if(p >= monoLen){
                    subs=subseq.substr(0,kseq_prime_len-(p-monoLen)-1);
                    pos=p+1-kseq_prime_len;
                }
                score= singleStringScore(score, &subs, pos, map);
                maximum_score=(score > maximum_score)?score:maximum_score;
            }
        }
    }
    else{
        for(int p=0; p < monoLen + kseq_prime_len -1; p++){
            std::string subs;
            double score=0;
            long pos;
            if( p < (monoLen-1)/2){
                if(p < kseq_prime_len-1){
                    subs=subseq.substr(kseq_prime_len-1-p,p+1); // for kmers longer than the 1/2 matrix length
                    pos=0;
                }
                else if(p >=kseq_prime_len-1){  //for kmers shorter than the 1/2marix length
                    subs=subseq;
                    pos=p-kseq_prime_len+1;
                }
                score= singleStringScore(score, &subs, pos, map);
                maximum_score=(score > maximum_score)?score:maximum_score;
            }
            else if(p >= (monoLen-1)/2  && (p - (monoLen-1)/2 <= kseq_prime_len -1)){
                std::string left_seq;
                long single_pos;
                if(kseq_prime_len <= monoLen && kseq_prime_len >= 0.5 * (monoLen-1)){
                    if(2*(p-(monoLen-1)/2 ) <=kseq_prime_len-1){
                        subs=subseq.substr(kseq_prime_len-2*(p+1-(monoLen-1)/2)+1,2*(p+1-(monoLen-1)/2)-1);
                        single_pos=p+1-subs.size();
                        pos=0;
                        if(p+1 >= kseq_prime_len){
                            pos=p+1-kseq_prime_len;
                            left_seq=subseq.substr(0,kseq_prime_len-2*(p-(monoLen-1)/2)-1); //if($kseq_prime_len-2*($    p-($mono_len/2 -1)) > 0);
                        }
                        else{
                            pos=0;
                            left_seq=subseq.substr(kseq_prime_len-(p+1),p+1-2*(p-(monoLen-1)/2)-1);
                        }
                    }
                    else{
                        subs=subseq.substr(0,2*(kseq_prime_len-(p+1-(monoLen-1)/2 ))+1);
                        pos= p-kseq_prime_len+1+subs.size();
                        single_pos=p+1-kseq_prime_len;
                        if(p >= monoLen){
                            std::string left_seq_pre=subseq.substr(2*(kseq_prime_len-(p+1-(monoLen-1)/2))+1);
                            left_seq=left_seq_pre.substr(0,kseq_prime_len-(p+1-monoLen)-subs.size());
                        }
                        else{
                            left_seq=subseq.substr(2*(kseq_prime_len-(p+1-(monoLen-1)/2))+1);
                        }
                    }
                }
                else if( kseq_prime_len > monoLen){
                    std::string left_seq_pre;
                    if(p < monoLen-1){
                        subs=subseq.substr(kseq_prime_len-2*(p+1-(monoLen-1)/2)+1,2*(p+1-(monoLen-1)/2)+1);
                        if(2*(p-(monoLen/2 -1)) < kseq_prime_len)
                            left_seq_pre=subseq.substr(0,kseq_prime_len-2*(p+1-(monoLen-1)/2)+1);
                        left_seq=left_seq_pre.substr(kseq_prime_len - (p+1));
                        pos=0;
                        single_pos=p+1-subs.size();
                        }
                    else if( p >= monoLen-1 && p <= kseq_prime_len-1){
                        subs=subseq.substr(kseq_prime_len-(p+1),monoLen);
                        pos=0;
                        single_pos=0;
                        }
                    else if(p > kseq_prime_len -1){
                        subs=subseq.substr(0,2*(kseq_prime_len-(p+1-(monoLen-1)/2 ))+1);
                        left_seq_pre=subseq.substr(2*(kseq_prime_len-(p+1-(monoLen-1)/2 ))+1);
                        left_seq=left_seq_pre.substr(0, left_seq_pre.size()-(p+1-monoLen));
                        pos=p-kseq_prime_len+1+subs.size();
                        single_pos=p+1-kseq_prime_len;
                    }
                }
                else if(kseq_prime_len <  0.5 * (monoLen-1)){
                    if(2*(p-(monoLen-1)/2 ) <=kseq_prime_len-1){
                        subs=subseq.substr(kseq_prime_len-2*(p+1-(monoLen-1)/2)+1,2*(p+1-(monoLen-1)/2)+1);
                        pos=p+1-kseq_prime_len;
                        left_seq=subseq.substr(0,kseq_prime_len-2*(p-(monoLen-1)/2)-1);
                        single_pos=p+1-subs.size();
                    }
                    else{
                        subs=subseq.substr(0,2*(kseq_prime_len-(p+1-(monoLen-1)/2))+1);
                        pos= p-kseq_prime_len+1+ subs.size();
                        left_seq=subseq.substr(2*(kseq_prime_len-(p+1-(monoLen-1)/2))+1);
                        single_pos=p+1-kseq_prime_len;
                    }
                }
                score=pairStringScore(score,&subs,single_pos,map);
                if(!left_seq.empty())
                    score= singleStringScore(score, &left_seq, pos, map);
                maximum_score=(score > maximum_score)?score:maximum_score;
            }
            else if(p >= (monoLen-1)/2 +kseq_prime_len -1){
                if(p <monoLen){
                    subs=subseq;
                    pos=p+1-kseq_prime_len;
                    }
                else if(p >= monoLen){
                    subs=subseq.substr(0,kseq_prime_len-(p-monoLen)-1);
                    pos=p+1-kseq_prime_len;
                }
                score= singleStringScore(score, &subs, pos, map);
                maximum_score=(score > maximum_score )?score:maximum_score;
            }
        }
    }
    return maximum_score;
}
//sub-routines for the dangling sequence
inline double singleStringScore(double x, std::string* s, long p, std::unordered_map<std::string, std::vector<double> > *m){
    for(int j=0; j<s->size(); j++){
        std::string letter=s->substr(j,1);
        if(letter == "-")
            continue;
        x+= (*m)[letter][p+j];
    }
    return x;
}
//sub-routines for the paired sequence
inline double pairStringScore(double x, std::string * s, long sp, std::unordered_map<std::string,std::vector<double> > *m){
    long subs_len=s->size();
    long diLen=(*m)["AA"].size();
    int half_subs_len;
    if(subs_len %2==0)
        half_subs_len=0.5* subs_len ;
    else
        half_subs_len=0.5* (subs_len-1) ;
    
    for (int j=0; j <= half_subs_len -1 ; j++){
        std::string letterA=s->substr(j,1);
        std::string letterB=s->substr(subs_len-j-1,1);
        std::string letter;

        if(letterA =="-" && letterB=="-")
            continue;
        else if(letterA =="-" ||letterB=="-"){
            long spos;
            if(letterB=="-"){
                spos=sp +j ;
                letter=letterA;
            }//#letterA is not gap so use the postion from the subs' beginning
            else if(letterA=="-"){
                spos=sp+(subs_len-1)-j; //#B is not gap, so from the end of subs
                letter=letterB;
            }
            x+=(*m)[letter][spos];
        }
        else{
            letter=letterA+letterB;
//            std::cout << *s << "\t" << letter  << "\t" << diLen-half_subs_len+j<< std::endl;
            x+=(*m)[letter][diLen-half_subs_len+j];
        }
    }
    return x;
}



