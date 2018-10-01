// =============================================================================
// CD-HIT-FP
//
// CD-HIT-FP is a very fast clustering program to cluster similar compounds 
// from a very large small molecule library.
//
// program written by
//                                      Weizhong Li
//                                      Email liwz@sdsc.edu
// Reference:
// Weizhong Li. A Fast Clustering Algorithm for Analyzing Highly Similar 
// Compounds of Very Large Libraries. J. Chem. Inf. Model. (2006) 46:1919-1923
// =============================================================================

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<time.h>

using namespace std;

#define MAX_REG_NAME 41 

int              max_fprint_no = 6000000;

double           threshold = 0.8;
int              long_bit = 1;  // 1: choose long Na as rep

int              fprint_no;
unsigned char *(*fprint);
int             *fprint_Na;
char          *(*fprint_name);
int             *fprint_sorted_idx;
int             *fprint_clstr_idx;
char            *fprint_iden;

int              clstr_no;
int             *clstr_rep;
int             *clstr_population;
int           *(*clstr_members);
int              fprint_seg = 0;
int           *(*fprint_seg_Na);


int len_fprintb= 988;
int len_fprint = 124;
char *type_fprint;

char fp_filename[256];
char clstr_filename[256];
char fp_out_filename[256];
int  fp_out_flag = 0;

#include "fp_class.h"

int  print_usage(char *arg);
int  fprint_write(ostream &out1, unsigned char *(*fp1),
                char *(*name1),  int no1); 

int main(int argc, char **argv) {
  int i, j, k;
  int i0, j0, k0;
  time_t current_time;
  if ((type_fprint=new char[80]) == NULL) bomb_error("memory");
  strcpy(type_fprint, "generic");

  if (argc <= 1) print_usage(argv[0]);

  for (i=1; i<argc; i++) {
    if (strcmp(argv[i], "-i") == 0)
      strcpy(fp_filename, argv[++i]);
    else if (strcmp(argv[i], "-o") == 0)
      strcpy(fp_out_filename, argv[++i]);
    else if (strcmp(argv[i], "-c") == 0)
      threshold = atof(argv[++i]);
    else if (strcmp(argv[i], "-l") == 0)
      long_bit = atoi(argv[++i]);
    else if (strcmp(argv[i], "-n") == 0)
      len_fprintb = atoi(argv[++i]);
    else if (strcmp(argv[i], "-t") == 0)
      strcpy(type_fprint, argv[++i]);
    else if (strcmp(argv[i], "-m") == 0) {
      max_fprint_no = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "-s") == 0) {
      fprint_seg = atoi(argv[++i]);
      if (fprint_seg<0 || fprint_seg==1) {print_usage(argv[0]); exit(1);}
    } 
    else if (strcmp(argv[i], "-f") == 0)
      fp_out_flag = atoi(argv[++i]);
    else
      strcpy(fp_filename, argv[i]);
  }



  len_fprint = len_fprintb/8;
  if (len_fprintb%8) len_fprint++; 

  strcpy(clstr_filename, fp_out_filename);
  strcat(clstr_filename, ".clstr");

  ofstream fclstr(clstr_filename);
  ofstream fout; if (fp_out_flag) fout.open(fp_out_filename);
  ifstream fpin(fp_filename);

  if ((fprint_name = new char *[max_fprint_no]) == NULL) 
    bomb_error("Memory");
  if ((fprint   = new unsigned char *[max_fprint_no]) == NULL)
     bomb_error("Memory");

  for (i=0; i<max_fprint_no; i++) {
    if ((fprint_name[i] = new char [MAX_REG_NAME]) == NULL)
      bomb_error("Memory");
    if ((fprint[i]   = new unsigned char [len_fprint]) == NULL)
      bomb_error("Memory");
    fprint_name[i][0] = 0;
    for (j=0; j<len_fprint; j++) fprint[i][j] = 0;
  }

  fprint_read(fpin, fprint, fprint_name, fprint_no, len_fprint, type_fprint);
  fpin.close();

  init_byte_2_bits();

  if ((fprint_Na = new int [fprint_no]) == NULL) bomb_error("Memory");
  if ((fprint_sorted_idx = new int [fprint_no]) == NULL) bomb_error("Memory");
  if ((fprint_clstr_idx = new int [fprint_no]) == NULL) bomb_error("Memory");
  if ((fprint_iden = new char [fprint_no]) == NULL) bomb_error("Memory");
  for (i=0; i<fprint_no; i++) fprint_Na[i] = 0;
  if (fprint_seg) {
    if ((fprint_seg_Na = new int *[fprint_no]) == NULL) bomb_error("Memory");
    for (i=0; i<fprint_no; i++) 
      if ((fprint_seg_Na[i] = new int [fprint_seg]) == NULL) 
        bomb_error("Memory");
    fprint_bit_count_seg(fprint_seg_Na, fprint_seg, 
                         fprint, fprint_no, len_fprint);
  }
  
  fprint_bit_count(fprint_Na, fprint, fprint_no, len_fprint);
  //fprint_sort(fprint_Na, fprint_sorted_idx, fprint_no);
  fprint_sort_uniq(fprint_Na, fprint, len_fprint, fprint_sorted_idx, fprint_no);
  if (long_bit) revert_array(fprint_sorted_idx, fprint_no);
  cout << fprint_no << " fprints read in" << endl;
//  fprint_write(cout, fprint, fprint_name, fprint_no);

  //for (i0=0; i0<fprint_no; i0++) {
  //  i = fprint_sorted_idx[i0];
  //  cout << fprint_name[i] << "\t" << fprint_Na[i] << endl;
  //}

  // begin clustering:
  unsigned char *this_fp;
  unsigned char *that_fp;
  char           this_iden;
  int            this_Na, that_Na;
  int            threshold_iden = (int) (threshold*100);

  clstr_no = 0;
  if ((clstr_rep = new int [fprint_no]) == NULL) bomb_error("Memory");

  for (i0=0; i0<fprint_no; i0++) {
    i = fprint_sorted_idx[i0];
    this_fp   = fprint[i];
    this_iden = 100;
    this_Na   = fprint_Na[i];

    // if >= 0, flag hold the clstr no this_fp similar to 
    int flag = -1;

    if (this_Na) { // ensure this_Na is not ZERO
      if (long_bit) { // that is longer than this
        for (j0=clstr_no-1; j0>=0; j0--) {
          j = clstr_rep[j0];
          that_fp = fprint[j];
          that_Na = fprint_Na[j];
          if (! that_Na) continue;
 
          // when following condition is ok, no need to calculate tanimoto
          if ( (this_Na*100/that_Na) < threshold_iden) break;
          if (fprint_seg) {
            this_iden = max_similarity_by_segs(fprint_seg, 
              fprint_seg_Na[j], fprint_seg_Na[i], that_Na, this_Na);
            if (this_iden < threshold_iden) continue;
          }
          this_iden = tanimoto_similarity(that_fp, this_fp, 
                                          that_Na, this_Na, len_fprint);
          if ( this_iden  >= threshold_iden ) {
            flag = fprint_clstr_idx[j]; break;
          }
        }
      }
      else { // not long_bit, then this is longer than that
        for (j0=clstr_no-1; j0>=0; j0--) {
          j = clstr_rep[j0];
          that_fp = fprint[j];
          that_Na = fprint_Na[j];
          if (! that_Na) continue;
    
          // when following condition is ok, no need to calculate tanimoto
          if ( (that_Na*100/this_Na) < threshold_iden) break;
          if (fprint_seg) {
            this_iden = max_similarity_by_segs(fprint_seg, 
              fprint_seg_Na[j], fprint_seg_Na[i], that_Na, this_Na);
            if (this_iden < threshold_iden) continue;
          }
          this_iden = tanimoto_similarity(that_fp, this_fp,
                                          that_Na, this_Na, len_fprint);
          if ( this_iden  >= threshold_iden ) {
            flag = fprint_clstr_idx[j]; break;
          }
        }
      }
    } // if this_Na

    if ( flag >=0 ) { // similar to an old 
      fprint_clstr_idx[i]   = -1 - flag;
      fprint_iden[i]        = this_iden;
    }
    else { // became new clstr
      fprint_clstr_idx[i]   = clstr_no;
      fprint_iden[i]        = 100;
      clstr_rep[clstr_no++] = i;
    }

    if ((i0+1) % 500 == 0 ) {
      cout << ".";
      if ((i0+1) % 5000 == 0) 
        cout << i0+1 << " fprints\t" << clstr_no << " clusters" << endl;
    } 
  } // for (i0=0; i0<fprint_no; i0++) {
  cout << i0 << " fprints\t" << clstr_no << " clusters" << endl;


  if ((clstr_population = new int [clstr_no]) == NULL) bomb_error("Memory");
  if ((clstr_members = new int *[clstr_no]) == NULL) bomb_error("Memory");
  for(i=0; i<clstr_no; i++) clstr_population[i] = 0;

  for(i=0; i<fprint_no; i++) {
    k = fprint_clstr_idx[i];
    if (k<0) k = -k - 1;
    clstr_population[k]++;
  }
  for(i=0; i<clstr_no; i++) {
    if ((clstr_members[i] = new int [clstr_population[i]]) == NULL) 
      bomb_error("Memory");
  }

  for(i=0; i<clstr_no; i++) clstr_population[i] = 0;
  //for(i=0; i<fprint_no; i++) {
  for (i0=0; i0<fprint_no; i0++) {
    i = fprint_sorted_idx[i0];
    k = fprint_clstr_idx[i];
    if (k<0) k = -k - 1;
    clstr_members[k][ clstr_population[k] ] = i;
    clstr_population[k]++;
  }

  // print out clstr info
  for(i=0; i<clstr_no; i++) {
    fclstr << ">Clstr " << i << " " << clstr_population[i] << endl;
    for (j=0; j<clstr_population[i]; j++) {
      k = clstr_members[i][j];
      k0 = fprint_iden[k];
      j0 = fprint_Na[k];
      fclstr << j << "\t" << fprint_name[k] << "\t" <<j0 << "\t";
      if (clstr_rep[i] == k) fclstr << "*" << endl;
      else                   fclstr << k0 << "%" << endl;
    }
  }
  fclstr.close();

  if (fp_out_flag) {
    ifstream fpin1(fp_filename);
    fprint_read_write(fpin1, fout, fprint_clstr_idx, len_fprint, type_fprint);
    fpin1.close();
    fout.close();
  }
  return 0;
}
// END main

int fprint_write(ostream &out1, unsigned char *(*fp1),
                 char *(*name1),  int no1) {
  int i0, i, j, k;
  unsigned char c1, c2, c3;

  for (i0=0; i0<no1; i0++) {
    i = fprint_sorted_idx[i0];
    out1 << name1[i] << "\t" << fprint_Na[i] << endl;
    for (j=0; j<len_fprint; j++) {
      c1 = fp1[i][j];
      c2 = c1 >> 4;
      c3 = c1 & 15;
      c2 = (c2 < 10) ? c2+'0' : c2-10+'a';
      c3 = (c3 < 10) ? c3+'0' : c3-10+'a';
      out1 << c2 << c3;
    }
    out1 << endl;
  }
  return 1;
}
// END fprint_write


int print_usage (char *arg) {
  cout << "Usage "<< arg << " [Options] \n";
  cout << "Options\n";
  cout << "\t-i input fingerprint file, required\n";
  cout << "\t-o output clstr file, required\n";
  cout << "\t-t type of input fingerprint, default generic\n";
  cout << "\t   can also be unity\n";
  cout << "\t-n length of a fingerprint, default 988, for unity fingerprint\n";
  cout << "\t-f 1 or 0, default 0,\n";
  cout << "\t   if set to 1, it write out fingerprint file for repsentatives\n";
  cout << "\t-c clustering threshold (Tanimoto coefficient), default 0.8\n";
  cout << "\t-l 1 or 0, default 1,\n";
  cout << "\t   choose fingerprints with most bit \"1\" as representative\n";
  cout << "\t   if set to 0, then choose fingerprints with most bit \"0\" \n";
  cout << "\t-m, maxium number of input compounds, default 6000000\n";
  cout << "\t    it should be >= real number of input compounds\n";
  cout << "\t-s, No. segment, default 0, \n";
  cout << "\t    if set to a positive number (2,3,4), it may increase the speed\n";
  cout << "\t    it depends on your input, you may need to try which is best\n";
  cout << "\t-h, print this help\n\n\n";
  cout << "\tif you find this program useful, please cite:\n";
  cout << "\tWeizhong Li. A Fast Clustering Algorithm for Analyzing Highly Similar Compounds\n";
  cout << "\tof Very Large Libraries. J. Chem. Inf. Model. (2006) 46:1919-1923\n\n\n";
  cout << "Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n\n";
  exit(1);
} // END print_usage



