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

int  init_byte_2_bits();
int  init_required_Nab(int *(*rNab), int high_bit, int threshold_iden);
void bomb_error(const char *message);
int  fprint_read(ifstream &in1, unsigned char *(*fp1),
                char *(*name1),  int &no1, int len_fp1, char *type_fprint);
int  fprint_read_write(ifstream &in1, ofstream &out1, int *fp_clstr_idxi,
                       int len_fp1, char *type_fprint);
int  fprint_read_generic(ifstream &in1, unsigned char *(*fp1),
                         char *(*name1),  int &no1, int len_fp1);
int  fprint_read_write_generic(ifstream &in1, ofstream &out1, 
                               int *fp_clstr_idxi, int len_fp1);
int  fprint_read_unity(ifstream &in1, unsigned char *(*fp1),
                       char *(*name1),  int &no1, int len_fp1);
int  fprint_read_write_unity(ifstream &in1, ofstream &out1, int *fp_clstr_idxi,
                             int len_fp1);
int  fprint_read_openbabel(ifstream &in1, unsigned char *(*fp1),
                         char *(*name1),  int &no1, int len_fp1);
int  fprint_read_write_openbabel(ifstream &in1, ofstream &out1, 
                               int *fp_clstr_idxi, int len_fp1);
int  fprint_bit_count(int *bit1, unsigned char *(*fp1), int no1, int len_fp1);
int  fprint_bit_count_seg(int *(*bit1), int seg1,
                          unsigned char *(*fp1), int no1, int len_fp1);
int  max_similarity_by_segs(int seg1, int *bit1, int *bit2, int Na, int Nb);
int  tanimoto_similarity(unsigned char *fp1, unsigned char *fp2,
                         int Na, int Nb, int len_fp1);
int  tanimoto_dissimilarity(unsigned char *fp1, unsigned char *fp2,
                            int Na, int Nb, int len_fp1);
int  fprint_sort(int *bit1, int *idx1, int no1);
int  fprint_sort0(int *bit1, int *idx1, int no1);
int  fprint_shuffling(int *idx1, int no1, int seed1);
int  quick_fprint_sort(int *a, int *idx,  int lo0, int hi0 );
int  quick_fprint_sort0(int *a, int *idx,  int lo0, int hi0 );
int  revert_array(int *a, int no1);

int  fprint_sort_uniq(int *bit1, unsigned char *(*fp1), int len_fp1, int *idx1, int no1);
int  quick_fprint_sort_uniq(int *a, unsigned char *(*fp1), int len_fp1, int *idx,  int lo0, int hi0 );
int  fpa_less_than_fpb(int Na, unsigned char *fp1, int Nb, unsigned char *fp2, int len_fp1);
int  fpa_big_than_fpb(int Na, unsigned char *fp1, int Nb, unsigned char *fp2, int len_fp1);

