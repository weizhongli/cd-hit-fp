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

using namespace std;

#define MAX_LINE_SIZE 10240

#include "fp_class.h"

int fprint_read(ifstream &in1, unsigned char **fp1, char **name1, int &no1, int len_fp1, char *type_fprint) {
  if (strcmp(type_fprint, "generic") == 0)      fprint_read_generic(in1, fp1, name1, no1, len_fp1);
  else if (strcmp(type_fprint, "unity") == 0)   fprint_read_unity(in1, fp1, name1, no1, len_fp1);
  else if (strcmp(type_fprint, "openbabel") == 0)   fprint_read_openbabel(in1, fp1, name1, no1, len_fp1);
  else bomb_error("unknown input format");
}
// END fprint_read

int fprint_read_generic(ifstream &in1, unsigned char **fp1, char **name1, int &no1, int len_fp1){
  char buffer1[MAX_LINE_SIZE];
  char *buffer2;
  int i, j, k;
  int len;
  int flag1;
  unsigned char c1;

  no1 = 0;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    len = strlen(buffer1);
    if (len == 0) continue;

    //get fingerprint name
    flag1 = 0;
    for (i=0; i<len; i++) {
      if (buffer1[i] == '\t') {
         buffer1[i] = 0; 
         flag1 = 1;
         break;
      }
    }
    if (! flag1 ) {cerr << buffer1; bomb_error(buffer1);}
    strcpy(name1[no1], buffer1);
    
    //get fingerprints and transform to 0 and 1 format
    buffer2=buffer1+i+1;
    len = strlen(buffer2);
    if (len != len_fp1*2) {
      cerr << "wrong length of fprint" << endl;
      bomb_error(buffer1);
    }
    for (i=0; i<len; i++) {
      c1 = buffer2[i];
      c1 = (c1 <= '9') ? c1 - '0' : c1 - 'a' + 10;
      j = i/2;
      fp1[no1][j] |= (i%2 == 0) ? c1 << 4 : c1;
    }

    no1++;
  }
  return 1;
}
// END fprint_read_generic


int fprint_read_unity(ifstream &in1, unsigned char **fp1, char **name1, int &no1, int len_fp1){
  char buffer1[MAX_LINE_SIZE];
  char buffer2[MAX_LINE_SIZE];
  int i, j, k;
  int len;
  int flag1;
  unsigned char c1;

  no1 = 0;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if (buffer1[0] == '#') continue;

    if ( isdigit(buffer1[0]) ) {
      len = strlen(buffer1);
      flag1 = 0;
      for (i=0; i<len; i++) {
        if (buffer1[i] == '\t') {
           buffer1[i] = 0; flag1 = 1;
        }
      }
      if (! flag1 ) bomb_error(buffer1);
      len = strlen(buffer1);
      for (i=0,j=0; i<len; i++) {
        if (buffer1[i] == ' ') {
          j++;
          if (j == 2) {
            strcpy(name1[no1], buffer1+i+1);
            break;
          }
        }
      }

      in1.getline(buffer2, MAX_LINE_SIZE-2, '\n');
      len = strlen(buffer2);
      if (len != len_fp1*2) {
        cerr << "wrong length of fprint" << endl;
        bomb_error(buffer1);
      }

      for (i=0; i<len; i++) {
        c1 = buffer2[i];
        c1 = (c1 <= '9') ? c1 - '0' : c1 - 'a' + 10;
        j = i/2;
        fp1[no1][j] |= (i%2 == 0) ? c1 << 4 : c1;
      }

      in1.getline(buffer2, MAX_LINE_SIZE-2, '\n');
      no1++;
    }
    else {
      continue;
    } 
  }
  return 1;
}
// END fprint_read_unity


int fprint_read_openbabel(ifstream &in1, unsigned char **fp1, char **name1, int &no1, int len_fp1){
  char buffer1[MAX_LINE_SIZE];
  char *buffer2;
  int i, j, k;
  int len;
  int flag1;
  unsigned char c1;

  no1 = -1;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    len = strlen(buffer1);
    if (len == 0) continue;

    if (buffer1[0] == '>') {
      //get fingerprint name
      for (i=0; i<len; i++) {
        if (isspace(buffer1[i])) {
           buffer1[i] = 0; 
           break;
        }
      }
      no1++;
      strcpy(name1[no1], buffer1+1);
      j=0; k=0;
    }
    else {
      for (i=0; i<len; i++) {
        c1=buffer1[i];
        if (isspace(c1)) continue;
        c1 = (c1 <= '9') ? c1 - '0' : c1 - 'a' + 10;
        j = k/2;
        fp1[no1][j] |= (k%2 == 0) ? c1 << 4 : c1;
        k++;
        if (k > len_fp1*2) {
          cerr << "wrong length of fprint" << endl;
          bomb_error(buffer1);
        }
      }
    }
  }
  no1++;
  return 1;
}
// END fprint_read_openbabel


//=================write cluster representatives=========================================
int fprint_read_write(ifstream &in1, ofstream &out1, int *fp_clstr_idx,
                      int len_fp1, char *type_fprint) {
  if      (strcmp(type_fprint, "generic") == 0)
    fprint_read_write_generic(in1, out1, fp_clstr_idx, len_fp1);
  else if (strcmp(type_fprint, "unity") == 0)
    fprint_read_write_unity(in1, out1, fp_clstr_idx, len_fp1);
  else if (strcmp(type_fprint, "openbabel") == 0)
    fprint_read_write_openbabel(in1, out1, fp_clstr_idx, len_fp1);
  else
    bomb_error("noknown input format");
}
// END fprint_read_write

int fprint_read_write_generic(ifstream &in1, ofstream &out1, int *fp_clstr_idx, int len_fp1) {
  char buffer1[MAX_LINE_SIZE];
  int i, j, k;

  int no1 = 0;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if (strlen(buffer1) == 0) continue;
    if (fp_clstr_idx[no1] >= 0) {
      out1 << buffer1 << "\n";
    }
    no1++;
  }
  return 1;
}
// END fprint_read_write_generic


int fprint_read_write_unity(ifstream &in1, ofstream &out1, int *fp_clstr_idx,
                      int len_fp1) {
  char buffer1[MAX_LINE_SIZE];
  char buffer2[MAX_LINE_SIZE];
  int i, j, k;
  int len;
  int flag1;
  unsigned char c1;


  int no1 = 0;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if (buffer1[0] == '#') { out1 << buffer1 << "\n";  continue; }

    if ( isdigit(buffer1[0]) ) {
      len = strlen(buffer1);
      flag1 = 0; for (i=0; i<len; i++) if (buffer1[i] == '\t') flag1 = 1;
      if (! flag1 ) bomb_error(buffer1);

      in1.getline(buffer2, MAX_LINE_SIZE-2, '\n');
      len = strlen(buffer2);
      if (len != len_fp1*2) {
        cerr << "wrong length of fprint" << endl;
        bomb_error(buffer1);
      }

      if (fp_clstr_idx[no1] >= 0) {
        out1 << endl;
        out1 << buffer1 << "\n";
        out1 << buffer2 << "\n";
      }

      in1.getline(buffer2, MAX_LINE_SIZE-2, '\n');
      if (fp_clstr_idx[no1] >= 0) out1 << buffer2 << "\n";

      no1++;
    }
    else {
      continue;
    } 
  }
  return 1;
}
// END fprint_read_write_unity


int fprint_read_write_openbabel(ifstream &in1, ofstream &out1, int *fp_clstr_idx, int len_fp1) {
  char buffer1[MAX_LINE_SIZE];
  int i, j, k;

  int no1 = -1;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if (buffer1[0] == '>')    no1++;
    if (fp_clstr_idx[no1] >= 0) {
      out1 << buffer1 << "\n";
    }
  }
  no1++;
  return 1;
}
// END fprint_read_write_openbabel

//==========================================================

int byte_2_bits[256];
int bits_Nab[256][256];

int init_byte_2_bits(){
  int i, j, k;
  unsigned char c1, c2, c3;

  for (i=0; i<256; i++) {
    c1 = 1;
    c3 = 0;
    for (k=0; k<8; k++) {
      if ( i & c1 ) c3++;
      c1 = c1 << 1;
    }
    byte_2_bits[i] = c3;
  }

  for (i=0; i<256; i++)
    for (j=0; j<256; j++)
      bits_Nab[i][j] = byte_2_bits[i & j];
  return 1;
}
// END init_byte_2_bits


int init_required_Nab(int *(*rNab), int high_bit, int threshold_iden) {
  int i, j, k, Nab1;

  if (threshold_iden==0) {
    for (i=0; i<=high_bit; i++)
      for (j=0; j<=high_bit; j++) 
        rNab[i][j] = 0;
     return 0;
  }
  else if (threshold_iden==100) {
    for (i=0; i<=high_bit; i++) 
      for (j=i; j<=high_bit; j++) 
        rNab[i][j] = rNab[j][i] = i;
    return 0;
  }
  else {
    for (i=0; i<=high_bit; i++) {
      for (j=i; j<=high_bit; j++) {
        Nab1 = threshold_iden * (i+j) / (100+threshold_iden) + 1;

        if (Nab1 > i) Nab1 = i; // if i == 0
        if (Nab1 == 0) {
          rNab[i][j] = rNab[j][i] = Nab1; continue;
        }
        for (k=Nab1; k>=0; k--) {
          if ( (100*k)/(i+j-k) >= threshold_iden ) Nab1 = k;
          else break;
        }
        rNab[i][j] = rNab[j][i] = Nab1;
      }
    }
    return 0;
  }

}// END init_required_Nab


int fprint_bit_count(int *bit1, unsigned char **fp1, int no1, int len_fp1) {
  for (int i=0; i<no1; i++) {
    int bb(0);
    for (int j=0; j<len_fp1; j++){
      bb += byte_2_bits[fp1[i][j]];
   }
    bit1[i] = bb;
  }
  return 1;
}
// END fprint_bit_count

//consecutive fprint_seg bytes were counted in different segments
int fprint_bit_count_seg(int **bit1, int seg1, unsigned char **fp1, int no1, int len_fp1) {
  int i, j, k, bb;
  for (i=0; i<no1; i++) {
    for (j=0; j<seg1; j++) bit1[i][j] = 0;
    for (j=0; j<len_fp1; j++) {
      k  = j % seg1;
      bit1[i][k] += byte_2_bits[fp1[i][j]];
    }
  }
  return 1;
}
// END fprint_bit_count_seg


int max_similarity_by_segs(int seg1, int *bit1, int *bit2, int Na, int Nb) {
  int Nab=0;
  for (int i=0; i<seg1; i++)
    Nab += (bit1[i] < bit2[i]) ? bit1[i] : bit2[i];
  return ( (100*Nab)/(Na+Nb-Nab));
}
// END max_similarity_by_segs

int tanimoto_similarity(unsigned char *fp1, unsigned char *fp2,
                        int Na, int Nb, int len_fp1) {
  int Nab=0;
  for (int i=0; i<len_fp1; i++){
    Nab += byte_2_bits[fp1[i] & fp2[i]];
  }
  return ( (100*Nab)/(Na+Nb-Nab));
}
// END tanimoto_similarity

int tanimoto_dissimilarity(unsigned char *fp1, unsigned char *fp2,
                        int Na, int Nb, int len_fp1) {
  int Nab=0;
  for (int i=0; i<len_fp1; i++) {
    Nab += byte_2_bits[fp1[i] & fp2[i]];
  }
  return (100 - (Nab*100)/(Na+Nb-Nab) );
}
// END tanimoto_dissimilarity


int fprint_shuffling(int *idx1, int no1, int seed1) {
  int i, j, k;
  int d1, d2, tmp;
  double dd;

  srand(seed1);

  for (i=0; i<no1; i++) idx1[i] = i;

  j = no1/2;
  for (i=0; i<j; i++) {
    k  = rand();
    dd = (double) k / RAND_MAX;
    d1 = (int) ((no1-1.0) * dd);

    k  = rand();
    dd = (double) k / RAND_MAX;
    d2 = (int) ((no1-1.0) * dd);

    tmp     = idx1[d1];
    idx1[d1] = idx1[d2];
    idx1[d2] = tmp;
  }
  return 0;
}
// END fprint_shuffling


//=================quick sort====================
// copied from fprint_sort, but it also sort fprints if Na is same
int  fprint_sort_uniq(int *bit1, unsigned char **fp1, int len_fp1, int *idx1, int no1){
  int i, j, k;
  for (i=0; i<no1; i++) idx1[i] = i;
  quick_fprint_sort_uniq(bit1, fp1, len_fp1, idx1, 0, no1-1);
  return 1;
}

// sort only idx
int quick_fprint_sort_uniq(int *a, unsigned char **fp1, int len_fp1, int *idx,  int lo0, int hi0 ) {
   int lo = lo0;
   int hi = hi0;
   int mid;
   unsigned char *mid_fp;
   int tmp;
   
   if ( hi0 > lo0) {
      mid    = a[ idx[(lo0 + hi0)/2] ];
      mid_fp = fp1[idx[(lo0 + hi0)/2] ];
      
      while( lo <= hi ) {
         while( (lo < hi0) && fpa_less_than_fpb(a[idx[lo]], fp1[idx[lo]], mid, mid_fp, len_fp1) ) lo++;
         while( (hi > lo0) && fpa_big_than_fpb(a[idx[hi]], fp1[idx[hi]], mid, mid_fp, len_fp1) ) hi--;
         if( lo <= hi ) {
            tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
            lo++; hi--;
         }
      }

      if( lo0 < hi ) quick_fprint_sort_uniq(a, fp1, len_fp1, idx, lo0, hi );
      if( lo < hi0 ) quick_fprint_sort_uniq(a, fp1, len_fp1, idx, lo, hi0 );
   }
   return 0;
}
// END quick_fprint_sort_uniq

int fpa_less_than_fpb(int Na, unsigned char *fp1, 
                      int Nb, unsigned char *fp2, int len_fp1) {
  int i, j, k;

  if (Na == Nb) {
    for (i=0; i<len_fp1; i++) {
      if (fp1[i] == fp2[i]) continue;
      return ((fp1[i] < fp2[i]) ? 1 : 0);
    }
    return 0;
  }
  else return ((Na < Nb) ? 1 : 0);
}
// END int  fpa_less_than_fpb


int fpa_big_than_fpb(int Na, unsigned char *fp1, 
                     int Nb, unsigned char *fp2, int len_fp1) {
  int i, j, k;

  if (Na == Nb) {
    for (i=0; i<len_fp1; i++) {
      if (fp1[i] == fp2[i]) continue;
      return ((fp1[i] > fp2[i]) ? 1 : 0);
    }
    return 0;
  }
  else return ((Na > Nb) ? 1 : 0);
}
// END int  fpa_big_than_fpb
//=================quick sort====================

// bit1 won't sort, but idx holding the index of sorted array
int fprint_sort(int *bit1, int *idx1, int no1) {
  int i, j, k;
  for (i=0; i<no1; i++) idx1[i] = i;
  quick_fprint_sort(bit1, idx1, 0, no1-1);
  return 1;
}
// END fprint_sort

// sort only idx, not a
int quick_fprint_sort(int *a, int *idx,  int lo0, int hi0 ) {
  int lo = lo0;
  int hi = hi0;
  int mid;
  int tmp;

  if ( hi0 > lo0) {
    mid = a[ idx[(lo0 + hi0)/2] ];
 
    while( lo <= hi ) {
      while( ( lo < hi0 ) && ( a[idx[lo]] < mid ) ) lo++;
      while( ( hi > lo0 ) && ( a[idx[hi]] > mid ) ) hi--;
      if( lo <= hi ) {
        tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
        lo++; hi--;
      }
    } // while

    if( lo0 < hi ) quick_fprint_sort(a, idx, lo0, hi );
    if( lo < hi0 ) quick_fprint_sort(a, idx, lo, hi0 );
  } // if ( hi0 > lo0)
  return 0;
}
// END quick_fprint_sort

// bit1 won't sort, but idx holding the index of sorted array
int fprint_sort0(int *bit1, int *idx1, int no1) {
  int i, j, k;
  for (i=0; i<no1; i++) idx1[i] = i;
  int *bit2;

  if ((bit2 = new int [no1]) == NULL) 
    bomb_error("memory");
  for (i=0; i<no1; i++) bit2[i] = bit1[i];

  quick_fprint_sort0(bit2, idx1, 0, no1-1);

  delete [] bit2;
  return 1;
}
// END fprint_sort0


// sort only idx, not a
int quick_fprint_sort0(int *a, int *idx,  int lo0, int hi0 ) {
  int lo = lo0;
  int hi = hi0;
  int mid;
  int tmp;

  if ( hi0 > lo0) {
    mid = a[ ( lo0 + hi0 ) / 2 ];
 
    while( lo <= hi ) {
      while( ( lo < hi0 ) && ( a[lo] < mid ) ) lo++;
      while( ( hi > lo0 ) && ( a[hi] > mid ) ) hi--;
      if( lo <= hi ) {
        tmp=a[lo];   a[lo]=a[hi];     a[hi]=tmp;
        tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
        lo++; hi--;
      }
    } // while

    if( lo0 < hi ) quick_fprint_sort0(a, idx, lo0, hi );
    if( lo < hi0 ) quick_fprint_sort0(a, idx, lo, hi0 );
  } // if ( hi0 > lo0)
  return 0;
}
// END quick_fprint_sort0


int revert_array(int *a, int no1) {
  int i, j, k;
  j = no1/2;

  for (i=0; i<j; i++) {
    k=a[i]; a[i]=a[no1-i-1]; a[no1-i-1]=k;
  }
  return 0;
}
// END revert_array

void bomb_error(const char *message) {
  cerr << "\nFatal Error\n";
  cerr << message << endl;
  cerr << "\nProgram halted !! \n\n";
  exit (1);
} // END void bomb_error
