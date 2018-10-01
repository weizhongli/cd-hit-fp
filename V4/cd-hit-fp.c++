// =============================================================================
// last version V3 apply large contiguous memory may not work
// so V4 apply for each fingerprint pointer
//
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
#include<omp.h>

using namespace std;

#define MAX_REG_NAME 41       //max length of fingerprint name

extern int byte_2_bits[256];

int ncpu=0;
int chunksize=10;
int len_fprintb= 988;      //fingerprint length in bits, default 988 for unity FPs
int len_fprint = 124;      //fingerprint bytes, default 124=(988+7)/8 for unity FPs
char type_fprint[80]="generic";      //default "generic", could be "unity"

char fp_filename[256];     //input file name
char clstr_filename[256];  //output cluster file name
char fp_out_filename[256]; //output representative fingerprint file name
int  fp_out_flag = 0;      //if output representative fingerprint

int      max_fprint_no = 6000000;   //maxium number of input compounds,  should be >= real number of input compounds
double   threshold = 0.8;           //Tanimoto coefficient, Tab=Nab/(Na+Nb-Nab): Nab=number of common 1s; Na,Nb=number of 1s in FPa and FPb
int      long_bit = 1;              //1: choose long Na as rep

int              fprint_no;            //real number of compounds
char          *(*fprint_name);         //fingerprint names, will contain max_fprint_no elements, each element has MAX_REG_NAME bytes
unsigned char *(*fprint);              //fingerprints, will contain max_fprint_no elements, each element has len_fprint bytes
//unsigned char	*clstrfprint;
unsigned char	**pclstrfprint;

int             *fprint_Na;
int             *fprint_sorted_idx;
int             *fprint_clstr_idx;     //cluster number to which the representative belongs, or members have value of (-1-clstrno), that is: 0cluster-1,1cluster-2,...
char            *fprint_iden;          //Tanimoto coefficient for each FP: 100 for representative, Tanimoto for cluster members

int              clstr_no;
int             *clstr_rep;
int             *clstr_population;
int           *(*clstr_members);
int              fprint_seg = 0;       //slice FP into 2,3,4 segments may increase speed
int           *(*fprint_seg_Na);



#include "fp_class.h"

int  print_usage(char *arg);
int  fprint_write(ostream &out1, unsigned char *(*fp1), char *(*name1),  int no1); 

int main(int argc, char **argv) {
	int i, k;
	int i0, j0, k0;
	time_t start_t,end_t;clock_t start, finish;

	if (argc <= 1) print_usage(argv[0]);
	for (i=0; i<argc; i++) cout<<argv[i]<<' ';cout<<endl;
	for (i=1; i<argc; i++) {
		if (strcmp(argv[i], "-i") == 0)       strcpy(fp_filename, argv[++i]);
		else if (strcmp(argv[i], "-o") == 0)  strcpy(fp_out_filename, argv[++i]);
		else if (strcmp(argv[i], "-c") == 0)  threshold = atof(argv[++i]);
		else if (strcmp(argv[i], "-C") == 0)  chunksize = atoi(argv[++i]);
		else if (strcmp(argv[i], "-l") == 0)  long_bit = atoi(argv[++i]);
		else if (strcmp(argv[i], "-n") == 0)  len_fprintb = atoi(argv[++i]);
		else if (strcmp(argv[i], "-t") == 0)  strcpy(type_fprint, argv[++i]);
		else if (strcmp(argv[i], "-m") == 0)  max_fprint_no = atoi(argv[++i]);
		else if (strcmp(argv[i], "-s") == 0)	{fprint_seg = atoi(argv[++i]);if(fprint_seg<0||fprint_seg==1){print_usage(argv[0]); exit(1);}} 
		else if (strcmp(argv[i], "-f") == 0)	fp_out_flag = atoi(argv[++i]);
		else if (strcmp(argv[i], "-T") == 0)	ncpu=atoi(argv[++i]);
		else												strcpy(fp_filename, argv[i]);
	}

#if defined (_OPENMP)
if(ncpu>0) omp_set_num_threads(ncpu);
#endif

	len_fprint = (len_fprintb+7)/8;

	if ((fprint_name = new char *[max_fprint_no]) == NULL)		bomb_error("Memory");
	for (i=0; i<max_fprint_no; i++) {
		if ((fprint_name[i] = new char [MAX_REG_NAME]) == NULL)	bomb_error("Memory");
		fprint_name[i][0] = 0;
	}
	
	if ((fprint = new unsigned char *[max_fprint_no]) == NULL)	bomb_error("Memory");
	for (i=0; i<max_fprint_no; i++) {
		if ((fprint[i] = new unsigned char [len_fprint]) == NULL)  bomb_error("Memory");
		memset(fprint[i],0,len_fprint);
	}
	
	if ((pclstrfprint = new unsigned char *[max_fprint_no]) == NULL)   bomb_error("Memory");
	for (i=0; i<max_fprint_no; i++) {
		if ((pclstrfprint[i] = new unsigned char [len_fprint]) == NULL)  bomb_error("Memory");
		memset(pclstrfprint[i],0,len_fprint);
	}
	
	init_byte_2_bits();   //init byte_2_bits[256] and bits_Nab[256][256], how many 1s in 0-255 and common 1s between them

	ifstream fpin(fp_filename);
	fprint_read(fpin, fprint, fprint_name, fprint_no, len_fprint, type_fprint);//read fprint_no fingerprints, store in fprint, fprint_name
	fpin.close();

	if ((fprint_Na = new int [fprint_no]) == NULL)         bomb_error("Memory");
	if ((fprint_sorted_idx = new int [fprint_no]) == NULL) bomb_error("Memory");
	if ((fprint_clstr_idx = new int [fprint_no]) == NULL)  bomb_error("Memory");
	if ((fprint_iden = new char [fprint_no]) == NULL)      bomb_error("Memory");
	for (i=0; i<fprint_no; i++) fprint_Na[i] = 0;
	if (fprint_seg) {
		if ((fprint_seg_Na = new int *[fprint_no]) == NULL) bomb_error("Memory");
		for (i=0; i<fprint_no; i++)
			if ((fprint_seg_Na[i] = new int [fprint_seg]) == NULL)   bomb_error("Memory");
		fprint_bit_count_seg(fprint_seg_Na, fprint_seg, fprint, fprint_no, len_fprint);//slice fingerprint to get number of 1s in segments
	}

	fprint_bit_count(fprint_Na, fprint, fprint_no, len_fprint);
	//fprint_sort(fprint_Na, fprint_sorted_idx, fprint_no);
	fprint_sort_uniq(fprint_Na, fprint, len_fprint, fprint_sorted_idx, fprint_no);
	if (long_bit) revert_array(fprint_sorted_idx, fprint_no);
	cout << fprint_no << " finger prints read and sorted" << endl;
//  fprint_write(cout, fprint, fprint_name, fprint_no);

	int *index_fp=new int[len_fprintb+1];
	int *index_clstr=new int[len_fprintb+1];
	for(i=0;i<=len_fprintb;i++) index_clstr[i]=0;

	int start_nbit=long_bit?len_fprintb:0;
	int Na,lastNa(-1);

	for (i0=0; i0<fprint_no; i0++) {
		i = fprint_sorted_idx[i0];
		Na=fprint_Na[i];
		if(Na != lastNa){
			if(long_bit)while(start_nbit>=Na) index_fp[start_nbit--]=i0;
			else while(start_nbit<=Na) index_fp[start_nbit++]=i0;
			lastNa=Na;
			//cout << fprint_name[i] << "\t" << fprint_Na[i] << '\t'<< i0<< endl;
		}
		//if(Na==0)cout<<fprint_name[i] << "\t" << fprint_Na[i] << '\t'<< i0<<'\t'<<i<< endl;
	}
	if(long_bit)while(start_nbit>=0) index_fp[start_nbit--]=i0;
	else while(start_nbit<=len_fprintb) index_fp[start_nbit++]=i0;

	strcpy(clstr_filename, fp_out_filename);
	strcat(clstr_filename, ".index");
	ofstream indextable(clstr_filename);
	indextable<<"#sorted fingerprints index table:\n";
	if(long_bit)	for(i=len_fprintb;i>=0;i--)indextable<<i<<'\t'<<index_fp[i]<<endl;
	else				for(i=0;i<=len_fprintb;i++)indextable<<i<<'\t'<<index_fp[i]<<endl;
	indextable<<i<<'\t'<<fprint_no<<endl;
	indextable.close();
	
	start_t=time(NULL);start=clock();

	// begin clustering:
	unsigned char *this_fp;
	int            this_iden;
	int            this_Na;
	int            threshold_iden = (int) (threshold*100);

	int j0end,clstr_Na_i=long_bit?len_fprintb:0;
	lastNa=-1;
	clstr_no = 0; 
	if ((clstr_rep = new int [fprint_no]) == NULL) bomb_error("Memory");
	for (i0=0; i0<fprint_no; i0++) {//clustering loop
		i = fprint_sorted_idx[i0];
		this_fp   = fprint[i];
		this_iden = 100;
		this_Na   = fprint_Na[i];
		int flag(-1);// if >= 0, flag hold the clstr no this_fp similar to 
		if (this_Na){ // ensure this_Na is not ZERO
			if (long_bit) { // that(cluster representative) is longer than this//calculate the start cluster number to compare
				j0end=this_Na*100/threshold_iden; j0end=j0end>len_fprintb?len_fprintb:j0end; j0end=index_clstr[j0end];
			}else{
				j0end=this_Na*threshold_iden/100; j0end=index_clstr[j0end];
			}

			#pragma omp parallel for schedule(static,10)
			for (j0=clstr_no-1; j0>=j0end; j0--) {
				if(flag>=j0){j0=-1;continue;}
				int j(clstr_rep[j0]);
				//unsigned char *that_fp(fprint[j]);
				unsigned char *that_fp(pclstrfprint[j0]);
				int that_Na(fprint_Na[j]);
				if (! that_Na) continue;

				// when following condition is true, no need to calculate tanimoto, this_Na/that_Na is the maximal Tanimoto value
				//if ( (this_Na*100/that_Na) < threshold_iden) j0=-1;// only this line is different between different long_bit
				int t_this_iden;
				if (fprint_seg) {
					t_this_iden = max_similarity_by_segs(fprint_seg, fprint_seg_Na[j], fprint_seg_Na[i], that_Na, this_Na);
					if (t_this_iden < threshold_iden) continue;
				}
				t_this_iden = tanimoto_similarity(that_fp, this_fp, that_Na, this_Na, len_fprint);
				if ( t_this_iden >= threshold_iden){
					if(j0>flag){
						#pragma omp critical
						{
							if(j0>flag) flag=j0;
						}
					}
				}
			}//end of for(j0)
		}//end of this_Na

		if ( flag >=0 ) { // similar to an old
			int j = clstr_rep[flag];
			unsigned char *that_fp = fprint[j];
			//unsigned char *that_fp = pfprint[j];
			int that_Na = fprint_Na[j];
			fprint_iden[i]= tanimoto_similarity(that_fp, this_fp, that_Na, this_Na, len_fprint);
			fprint_clstr_idx[i]   = -1 - flag;// i-th fingerprint similar to n-th cluster, value is (-1-n)
		}else{ // became a n-e-w clstr
			fprint_clstr_idx[i]   = clstr_no;
			fprint_iden[i]        = 100;

			//update cluster index table
			if(this_Na != lastNa){
				if(long_bit)while(clstr_Na_i>=this_Na) index_clstr[clstr_Na_i--]=clstr_no;
				else while(clstr_Na_i<=this_Na) index_clstr[clstr_Na_i++]=clstr_no;
				lastNa=this_Na;
			}
			memcpy(pclstrfprint[clstr_no],fprint[i],len_fprint);
			clstr_rep[clstr_no++] = i;
		}

		if ((i0+1) % 500 == 0 ) {
			cout << ".";
			if ((i0+1) % 5000 == 0) 
			cout << i0+1 << " fprints\t" << clstr_no << " clusters" << endl;
		}
	}//end of clustering loop
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
	for (i0=0; i0<fprint_no; i0++) {
		i = fprint_sorted_idx[i0];
		k = fprint_clstr_idx[i];
		if (k<0) k = -k - 1;
		clstr_members[k][ clstr_population[k] ] = i;
		clstr_population[k]++;
	}

	end_t=time(NULL);finish=clock();
	cout<<"Time used(seconds):\t"<<difftime(end_t,start_t)<<'\t'<<finish-start<<endl;

	// print out clstr info
	strcpy(clstr_filename, fp_out_filename);
	strcat(clstr_filename, ".clstr");
	ofstream fclstr(clstr_filename);
	for(i=0; i<clstr_no; i++) {
		fclstr << ">Clstr " << i << " " << clstr_population[i] << endl;
		for (int j=0; j<clstr_population[i]; j++) {
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
		ofstream fout(fp_out_filename);
		ifstream fpin1(fp_filename);
		fprint_read_write(fpin1, fout, fprint_clstr_idx, len_fprint, type_fprint);
		fpin1.close();
		fout.close();
	}
	
	strcpy(clstr_filename, fp_out_filename);
	strcat(clstr_filename, ".index");
	indextable.open(clstr_filename,fstream::app);
	indextable<<"#cluster index table:\n";
	if(long_bit)for(i=len_fprintb;i>=0;i--)indextable<<i<<'\t'<<index_clstr[i]<<endl;
	else			for(i=0;i<=len_fprintb;i++)indextable<<i<<'\t'<<index_clstr[i]<<endl;
	indextable<<i<<'\t'<<clstr_no<<endl;
	indextable.close();




	for(i=0; i<clstr_no; i++) {
		delete[] clstr_members[i];
	}
	delete[] clstr_members;	delete[] clstr_population;
	delete[] clstr_rep;
	delete[] index_fp;	delete[] index_clstr;
	delete[] fprint_Na;
	delete[] fprint_sorted_idx;
	delete[] fprint_clstr_idx;
	delete[] fprint_iden;
	if (fprint_seg) {
		for(i=0; i<fprint_no; i++)
			delete[] fprint_seg_Na[i];
		delete[] fprint_seg_Na;
	}
	for (i=0; i<max_fprint_no; i++) {
		delete[] fprint_name[i];
		delete[] fprint[i];
		delete[] pclstrfprint[i];
	}
	delete[] fprint_name;
	delete[] fprint;
	delete[] pclstrfprint;

	return 0;
}
// END main

int fprint_write(ostream &out1, unsigned char *(*fp1), char *(*name1), int no1) {
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

