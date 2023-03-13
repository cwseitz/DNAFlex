/* Define functions */
#define min(a, b) (((a) + (b) - fabs((a) - (b))) * 0.5)
/* End Define */


/* Headers from DNA Flex */
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "dna_flex.h"	// Functions used in DNA Flex
#include "umbrella.h"	// Functions needed for Cartesian reconstruction
#include "energy.h"	// Functions for calculating int energies

/*
#include "schna_ar_p.h"
#include <time.h>
/* End Headers */


/*--------------------------Beginning of DNAFlex-------------------------*/

int main(int argc, char *argv[])
{

int mult_folders = 0;	//1: input number of folder when executing ./MC_Chromatin 1; 2: no input when executing

char nucl_file[200];
sscanf(argv[1], "%s", nucl_file);	// get nucleosome positions
printf("%s\n",nucl_file);

int no_nucl;
get_num_nucl(&no_nucl,nucl_file);
int *ins_nucl;
ins_nucl = malloc(no_nucl*sizeof(int));
get_nucl_pos(no_nucl,ins_nucl,nucl_file);
printf("no_nucl %d\n",no_nucl);


char seq_file[200];
sscanf(argv[2], "%s", seq_file);	// get sequence
printf("%s\n",seq_file);


int buf_num_it;
buf_num_it = atoi (argv[3]);	// get number of structures to generate
//printf("%d\n",buf_num_it);


char conf_file[200];
sscanf(argv[4], "%s", conf_file);	// get configuration file for linker DNA
printf("%s\n",conf_file);

char out_folder[200];
sscanf(argv[5], "%s", out_folder);	// get output folder
printf("%s\n",out_folder);


char out[] = "/usr/people/jwalther/Programs/Chromatin/src_te";	// just change this in case mult_folders = 1 (it runs from 1 to 10 in subfolders)
//char out_folder[strlen(out)];	 //char out_folder[strlen(out)+2];	when mult_folders = 1 to attach another output folder

if(mult_folders == 1)
{
int fo_ind;
sscanf(argv[1], "%d", &fo_ind);	// for outputing in different folders
//sprintf(out_folder,"%s/%d",out,fo_ind);
//printf("%s\n",out_folder);
}
else{
					
char stif_file[] = "stif_chrom_nfr_half_flex_200.dat";

/* Initialize all variables at beginning of code!! */ 

int nmc,i0,ibuff2,juega,last,psi;
int iflag; 	//it should be a bool, but I dont know how to assign true and false to a bool in c
char **seq1, **ist; //seq1[1000][10], ist[1000][10];
double ***stif, ***stiff, **geom; // stif[1000][6][6],stiff[1000][6][6],geom[1000][6];
int ord[6];
double **xconfi, **xconf0, **xconf; // xconfi[6][1000],xconf0[6][1000],xconf[6][1000];
double **xconft, *ener; // xconft[6][1000],ener[1000]
double sc[6],dd[6];
double **dx; //dx[6][1000];
char *names[6];   				// No maximum size of char in C
char line;
double boltz=0.00198717;			// in kcal/mol/K
double enex;
int itot;
int i,j,k;	// Iteration parameters
int it_mc;	// Test parameter for # MC moves
char **output;  // output[1000][100];			// Writing helical coordinates in this string
int it_prog,phi;
double ind_energy;
double *ind_ener;
float sum_energy;
int selection;
double ***stif_tetra, ***stif_di, **geom_tetra, **geom_di; //stif_tetra[1000][6][6],stif_di[1000][6][6],geom_tetra[1000][6],geom_di[1000][6];
char **seq1_tetra, **seq1_di;  //seq1_tetra[1000][10], seq1_di[1000][10];
int it_snapshots, it_mcsteps, alpha, it_mod;
char **out_part;
double endenergy;

double **xyz_float;
int msec,msec1,msec2,msec3;
double **center_bp;
double r_bp;
char **pdb_data;
char **tmp_at;
char *baset;
char *strand;
int *base_num;
double **tmp_xyz;
double *q_fl;
int num_bp;
double **d_bp;
double **d_bp_fl;
double **lj_bp, **lj_bp_fl;
double **mid_bp;
double **dir_bp;
double **dir_ph;
double **pos_ph, **buf_pos_ph;
double **dh_ph;
double **dh_ph_fl;
double q_ph;
double **d_ph1, **d_ph2;
double **d_ph_fl;
double tot_dh_lj, tot_lj, tot_dh;
int k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13;
char **pdb_tot;
int tot_num_at;
int num_fl;	// Number of floating patricles
double shift;
double shift_in;
int acc;
int acc_ph;
double bp_excl_vol;
double bp_fl_excl_vol;
double lj_k,lj_cut_up,lj_cut_low,lj_a_bp,lj_a_fl;
double eps,kappa,f_pi_eps0,dh_cut_low,dh_cut_up;
double **xconf_buf;
double tot_en,prev_tot_en;
double d_ph_midbp;
double **xyz_float_buf;
int count_no_excl;


double ***orient_old, ***orient_new, **pos_new, ***orient_store, **pos_store, ***rot_store, ***rot_new;
double **pos_old;

int no_rot_pars = 3;
int vol_olap_start = 0;

char **ist_nucl, **ist_all, **ist_mc;
char **seq_tetra, **seq_tetra_gen;
double **geom_tetra_gen, ***stif_tetra_gen, **xconf_all;
double **xconf_mc, **xconf0_mc, ***stiff_mc, **xconfi_mc, **xconft_mc, **xconf_nucl;
double ** nucl_pos, *center_nucl, **tr_mat, *buf_vec;
double nucl_nucl_excl_vol, nucl_bp_excl_vol, nucl_fl_excl_vol;
double **d_nucl, **d_nucl_bp, **d_nucl_fl, **d_nucl_ph;
char **table_bp, **table_bp_triad, **table_nucl, **table_pos;
char filename[300];
int err_count, metr_on, acc_vol, acc_mc, start_sim;
acc_vol=0; acc_mc=0;
int len_seq;	// Variable to determine dynamic size of all arrays
int seq_pcs;	// Length of sequence pieces to save (i.e. 4 for tetramer seq)
int no_helpars;
int len_out_line;
int len_bp_nucl,len_seq_nucl;
int len_all,bp_all;
int itot_mc, real_itt, buf_real_itt;
//int no_nucl = sizeof(ins_nucl)/sizeof(ins_nucl[0]);
double buf_tot_en = 1000000;
double buf_enertot = 1000000;
double tmax,t;
int do_sim_ann, no_mc_move, get_itt;
double factor;
int rnd_start;
int *it_dna_mc, *it_acc_vol, *it_acc_mc, sel;
double ***linker_orient, **linker_pos, ***nucl_dna_orient, **nucl_dna_pos;
int *link_nucl_index, *ph_link_nucl_index;
double q_link, q_nucl_dna, *q_nucl, *arr_q_ph, cut_d_bp_nucl, cut_d_ph_nucl;
double **lj_nucl, **lj_nucl_bp, **lj_nucl_fl;
double lj_k_nucl, lj_k_nucl_bp, lj_k_nucl_fl, lj_a_nucl, lj_a_bp_nucl, lj_a_nucl_fl, lj_nucl_cut_low, lj_nucl_cut_up;
int m_nucl_bp,m_nucl_ph;
double e_nucl_lj, e_nucl_dh, e_nucl;
int ***in_q_ph_ph, **in_q_nucl_ph, **in_q_ph_ph_MIN;
double dh_nucl_cut_low, dh_nucl_cut_up, **dh_nucl, **dh_nucl_ph, **dh_nucl_fl;
int incr, red_bp_all, red_ph;
double e_lj_bp_bpfl, e_dh_ph_phfl;
double ***store_R_T_mst, ***store_R_T_i, ***R_T_mst, ***R_T_i, ***buf_R_T_mst, ***buf_R_T_i;
int it_en, want_to_save_en, int_en_save, detect_indep, it_each_en;
double **save_en, sh_min, **save_each_en, sed_co, rg;
int len_seg, *ind_linker, nucl_ind, min_ph;

int_en_save = 1000;		// Interval to save energy value
len_bp_nucl = 146;
seq_pcs = 10;
no_helpars = 6;
len_out_line = 100;

num_fl = 0;			// Number of floating particles

/* ------------- Accuracy of model ------------- */
acc = 7;			// bp-steps for calculating distance (usually 10); minimum I can go is 7 because of overlap in nucleosome
acc_ph = 2*acc;			// phosphate steps for calculating distance (usually 20); 2*acc to guarantee that every acc bp the phosphate distance is evaluated; should be even to guarantee coverage of both strands


get_seq_len(&len_seq,&itot,seq_file);

printf("len_seq %d \n",len_seq);
printf("no_nucl %d\n", no_nucl);

itot_mc = itot + no_nucl;
printf("itot_mc %d\n",itot_mc);

len_all = no_nucl*len_bp_nucl+itot_mc;
printf("len_all %d\n",len_all);

bp_all = len_all + 1;

red_bp_all = (int) (double)(bp_all-1)/(double)acc;
printf("red_bp_all %d\n",red_bp_all);
red_ph = (int) (double)(2*bp_all-2) / (double)acc_ph;
printf("red_ph %d\n",red_ph);


int alc_space;
if(len_seq < 256){alc_space = 256;}else{alc_space = len_seq;}

/* Allocate memory for arrays */
seq1 = char_aloc2d(alc_space,seq_pcs); 
ist  = char_aloc2d(len_seq,seq_pcs); 
stif = dbl_aloc3d(alc_space,no_helpars,no_helpars);
stiff = dbl_aloc3d(len_seq,no_helpars,no_helpars);
geom = dbl_aloc2d(alc_space,no_helpars);
xconfi = dbl_aloc2d(no_helpars,len_seq);
xconf0 = dbl_aloc2d(no_helpars,len_seq);
xconf = dbl_aloc2d(no_helpars,len_seq);
xconft = dbl_aloc2d(no_helpars,len_seq);
ener = malloc(itot_mc*sizeof(double));
dx = dbl_aloc2d(no_helpars,itot_mc);
output  = char_aloc2d(len_all+5,len_out_line);

seq1_tetra = char_aloc2d(256,seq_pcs);
seq1_di = char_aloc2d(16,seq_pcs);
stif_tetra = dbl_aloc3d(256,no_helpars,no_helpars);
stif_di = dbl_aloc3d(16,no_helpars,no_helpars);
geom_tetra = dbl_aloc2d(256,no_helpars);
geom_di = dbl_aloc2d(16,no_helpars);

out_part  = char_aloc2d(len_seq,len_out_line);
xyz_float = dbl_aloc2d(num_fl,no_helpars);
q_fl = malloc(len_seq*sizeof(double));
center_bp = dbl_aloc2d(bp_all,no_helpars);
/*pdb_data = char_aloc2d(2*len_seq,len_out_line);
tmp_at = char_aloc2d(2*len_seq,no_helpars);
baset = malloc(2*len_seq*sizeof(char));
strand = malloc(2*len_seq*sizeof(char));
base_num = malloc(2*len_seq*sizeof(int));
tmp_xyz = dbl_aloc2d(2*len_seq,no_helpars);*/ 
d_bp = dbl_aloc2d(red_bp_all,red_bp_all); 
d_bp_fl = dbl_aloc2d(bp_all,bp_all);	//printf("Test\n");	
lj_bp = dbl_aloc2d(red_bp_all,red_bp_all);	//printf("Test\n");	
lj_bp_fl = dbl_aloc2d(bp_all,bp_all);		
//printf("Test1\n");
mid_bp = dbl_aloc2d(bp_all,no_helpars);
dir_bp = dbl_aloc2d(bp_all,no_helpars);
dir_ph = dbl_aloc2d(bp_all,no_helpars);
pos_ph = dbl_aloc2d(2*bp_all,no_helpars);
buf_pos_ph = dbl_aloc2d(2*bp_all,no_helpars);
d_ph1 = dbl_aloc2d(2*red_ph,2*red_ph);
d_ph2 = dbl_aloc2d(red_ph,red_ph);
d_ph_fl = dbl_aloc2d(red_ph,bp_all);	
dh_ph = dbl_aloc2d(2*red_ph,2*red_ph); //printf("Test2\n");
dh_ph_fl = dbl_aloc2d(2*red_ph,bp_all);
//printf("Test3\n");
xconf_buf = dbl_aloc2d(no_helpars,itot_mc);
xyz_float_buf = dbl_aloc2d(num_fl,no_helpars);

pdb_tot = char_aloc2d(20*len_seq,len_out_line);
//printf("Test4\n");
orient_old = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
orient_new = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
orient_store = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
pos_old = dbl_aloc2d(len_all,no_rot_pars);
pos_new = dbl_aloc2d(len_all,no_rot_pars);
pos_store = dbl_aloc2d(len_all,no_rot_pars);
rot_new = dbl_aloc3d(2*len_all,no_rot_pars,no_rot_pars);
store_R_T_mst = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
store_R_T_i = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
R_T_mst = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
R_T_i = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
buf_R_T_mst = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
buf_R_T_i = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
//printf("Test5\n");
xconf_nucl = dbl_aloc2d(no_helpars,len_bp_nucl);
ist_nucl  = char_aloc2d(len_bp_nucl+1,seq_pcs);
ist_all  = char_aloc2d(len_all,seq_pcs);
xconf_all  = dbl_aloc2d(no_helpars,len_all);
seq_tetra_gen = char_aloc2d(256,seq_pcs);
stif_tetra_gen = dbl_aloc3d(256,no_helpars,no_helpars);
geom_tetra_gen = dbl_aloc2d(256,no_helpars);
//printf("Test6\n");

xconf_mc = dbl_aloc2d(no_helpars,itot_mc);
xconf0_mc = dbl_aloc2d(no_helpars,itot_mc);
stiff_mc = dbl_aloc3d(itot_mc,no_helpars,no_helpars);
xconfi_mc = dbl_aloc2d(no_helpars,itot_mc);
xconft_mc = dbl_aloc2d(no_helpars,itot_mc);
ist_mc  = char_aloc2d(itot_mc,seq_pcs);
//printf("Test7\n");
nucl_pos = dbl_aloc2d(no_nucl,3);
center_nucl = malloc(3*sizeof(double));
tr_mat = dbl_aloc2d(3,3);
buf_vec = malloc(3*sizeof(double));
d_nucl = dbl_aloc2d(no_nucl,no_nucl);
d_nucl_bp = dbl_aloc2d(no_nucl,red_bp_all);
d_nucl_fl = dbl_aloc2d(no_nucl,bp_all);
d_nucl_ph = dbl_aloc2d(no_nucl,2*red_ph);
//printf("Test8\n");
table_bp  = char_aloc2d(bp_all+5,len_out_line);
table_pos  = char_aloc2d(bp_all+no_nucl+5,len_out_line);
table_bp_triad  = char_aloc2d(4*bp_all+5,len_out_line);
table_nucl  = char_aloc2d(no_nucl+5,len_out_line);
//printf("Test9\n");
linker_orient = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
linker_pos = dbl_aloc2d(len_all,no_rot_pars);
nucl_dna_orient = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
nucl_dna_pos = dbl_aloc2d(len_all,no_rot_pars);
link_nucl_index = malloc((2*no_nucl+2)*sizeof(int));
ph_link_nucl_index = malloc((2*no_nucl+2)*sizeof(int));
arr_q_ph = malloc((2*len_all)*sizeof(double));
//printf("Test10\n");
lj_nucl = dbl_aloc2d(no_nucl,no_nucl);
lj_nucl_bp =  dbl_aloc2d(no_nucl,red_bp_all);
lj_nucl_fl =  dbl_aloc2d(no_nucl,bp_all);
//printf("Test11\n");
in_q_ph_ph = int_aloc3d(2*red_ph,2*red_ph,2); //printf("Test12\n");
in_q_ph_ph_MIN = int_aloc2d(red_ph,red_ph);
in_q_nucl_ph = int_aloc2d(no_nucl,2*red_ph);
q_nucl = malloc(no_nucl*sizeof(double));
ind_linker = malloc(2*(no_nucl+1)*sizeof(int));

dh_nucl = dbl_aloc2d(no_nucl,no_nucl);
dh_nucl_ph =  dbl_aloc2d(no_nucl,2*red_ph);
dh_nucl_fl =  dbl_aloc2d(no_nucl,bp_all);
//printf("Test13\n");
/* Initialize arrays and files for huge tables */

char **table_heli_shif;
char **table_heli_slid;
char **table_heli_rise;
char **table_heli_tilt;
char **table_heli_roll;
char **table_heli_twis;

table_heli_shif  = char_aloc2d(bp_all,len_out_line);
table_heli_slid  = char_aloc2d(bp_all,len_out_line);
table_heli_rise  = char_aloc2d(bp_all,len_out_line);
table_heli_tilt  = char_aloc2d(bp_all,len_out_line);
table_heli_roll  = char_aloc2d(bp_all,len_out_line);
table_heli_twis  = char_aloc2d(bp_all,len_out_line);

sprintf(filename,"%s/output_tables_helpar/table_all_shif.dat",out_folder);
FILE *file_shif;
file_shif=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_slid.dat",out_folder);
FILE *file_slid;
file_slid=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_rise.dat",out_folder);
FILE *file_rise;
file_rise=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_tilt.dat",out_folder);
FILE *file_tilt;
file_tilt=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_roll.dat",out_folder);
FILE *file_roll;
file_roll=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_twis.dat",out_folder);
FILE *file_twis;
file_twis=fopen(filename, "w");

/* ------- End Initialze for tables and files --------- */


srand (getpid());		// initialize random seed (the random seed is initialized to a value representing the current time (calling time) to generate a different value every time the program is run.)


clock_t start_all = clock(), diff_all; //Timer

/* icuan = 1;
icnt = 1; */


/* This is the part which should be read from a txt file (see fortran code) */

nmc=100000;
double agr=1.0;
int ibuff=10000;
ibuff2=100;
double temp=298.0,fact=0.3;
last=1000000;
sc[0]=5.148;
sc[1]=3.558;
sc[2]=2.052;
sc[3]=29.472;
sc[4]=34.266;
sc[5]=53.862;
names[0]="shif";
names[1]="slid";
names[2]="rise";
names[3]="tilt";
names[4]="roll";
names[5]="twis";

/* ----------------------------------End of part--------------------------------------- */


/* ------------------------------- Looping and selection variables ------------------------------- */


selection = 3;			// 0 = individual, 1 = dimer, 2 = tetramer, 3 = tetramer with dimer end


it_mcsteps = 1;		// number of MC steps in floating process before looking at position of particle
it_snapshots = 1;		// number of pdb's generated during the floating process
it_prog = 1;			// number of executing all the MC and reconstruction part (THE WHOLE PROGRAM)
it_mc = 1000*buf_num_it*it_mcsteps;	// number of MC steps to bring structure to minimum energy before particle inserted
				/* 100000 -> 56mer
				   10000000 -> 1000mer
				   2000000 -> 450mer
				   5000000 -> 750mer
				*/
metr_on = 2;			// 1 = DNA gets sampled with metropolis, 2 = all MC moves get accepted

do_sim_ann = 2;			// Do simmulated annealing (do_sim_ann = 1) instead of Metropolis for whole structure
tmax = 5000;			// Maximum temperature simmulated annealing
factor = 0.95;			// Factor to decrease temperature

rnd_start = 2;			// Use random starting structure (1) or use xconfi as starting structure (2) or step-by-step start (3)
min_ph = 1;			// Calculate DH with (1) minimum distance of phosphate pair of a bp; (2) with all combinations of the phosphates of each bp (= 4 interactions per bp/bp DH)

incr = 2;			// Increase acc for each time building a new structure: (1) = yes
want_to_save_en = 1;		// Save energy every int_en_save accepted moves: (1) = yes
detect_indep = 0;		// Save all energy values to detect after what # of structures the structures get independent


/* ------------------------------- END Looping and selection variables ------------------------------- */


/* ------------------------------- Parameters for dynamical part ------------------------------------- */


shift_in = 600.0;		// Volume occupied by floating particles
shift = 5.0;			// Range of Shift of floating particle per MC move

/* ---- Parameters for excluded volume ---- */
bp_excl_vol = 20;		// minimum distance between two bp centers (in Angstrom)
bp_fl_excl_vol = 10;		// minimum distance between bp center and floating particle (in Angstrom)

d_ph_midbp = 10.25;		// Distance of Phosphate of helical axis



/* ---- Parameters for LJ ---- */	
lj_k = 1.0*0.0257; 		// in eV (0.0257 eV = 1kT, value 1*0.0257  adapted from Arya 2014)
lj_cut_up = 200;		// 1.5*bp_excl_vol (= 30)
lj_cut_low = 0.1;

sh_min = 1.55;

lj_a_bp = bp_excl_vol/sh_min;	// = 10A as excluded volume (radius of excluded volume around center_bp pos)
lj_a_fl = bp_fl_excl_vol;	// 



/* ---- Parameters for DH ---- */
eps = 80.0;
kappa = 0.103; 			// unit: 1/A	taken from Arya 2014 paper (usually 0.033 for 10mM salt) Mueller 0.103 for 100mM salt
f_pi_eps0 = 1.1126*pow(10,-10);
q_ph = -1.0;

dh_cut_low = 0.1;
dh_cut_up = 200;			// bp_excl_vol (= 20.0)


/* ---- Define center of nucleosomes ---- */
center_nucl[0] = 48.84575510;	// Derived from the mean of bp positions of 1kx5 structure
center_nucl[1] = -0.01046939;
center_nucl[2] = 27.70510204;

nucl_nucl_excl_vol = 55;	// 2*5.5nm (since nucleosome is 11nm in diameter)	(DNA in nucleosome still considered seperately)
nucl_bp_excl_vol = 37.5;		// minimum distance between nucleosome and center of bp (55A+10A)
nucl_fl_excl_vol = 55;		// minimum distance between nucleosome and floating particle (radius of nucleosome = 55A)

/* ---- Partial charges phosphates ---- */
q_link = -1.0;			// linker DNA
q_nucl_dna = -1.0;		// Nucleosomal DNA


/* ---- Parameters for nucleosome ---- */
cut_d_bp_nucl = 10;		// cutoff # bp before and after nucleosome not to account for distance evaluation
cut_d_ph_nucl = 2*cut_d_bp_nucl;// cutoff # ph before and after nucleosome not to account for distance evaluation

int mult_nucl_en = 1;
if(min_ph == 1){
for(i=0;i<no_nucl;i++){q_nucl[i] = mult_nucl_en*94/acc_ph;}	// charge histone core (2*146*q_nucl_dna + q_nucl = -198 = total charge of nucl core particle)
		}
else{
for(i=0;i<no_nucl;i++){q_nucl[i] = mult_nucl_en*94/(acc_ph/2);}	// charge histone core (2*146*q_nucl_dna + q_nucl = -198 = total charge of nucl core particle)
     }

/* - LJ nucl - */
lj_nucl_cut_low = 0.1;
lj_nucl_cut_up = 200;

lj_k_nucl = 1*lj_k;		// random
lj_k_nucl_bp = lj_k;		// random
lj_k_nucl_fl = lj_k;
//lj_a_nucl = nucl_nucl_excl_vol+55;
lj_a_nucl = 110/sh_min;
lj_a_bp_nucl = nucl_bp_excl_vol/sh_min;
lj_a_nucl_fl = nucl_fl_excl_vol/sh_min;

/* - DH nucl - */
dh_nucl_cut_low = 0.1;
dh_nucl_cut_up = 500;


/* ------------------------------- END Parameters for dynamical part ----------------------------------- */

int num_red_mc = it_mc/it_mcsteps;
int no_save_en = num_red_mc/int_en_save + 1;

it_dna_mc = malloc(num_red_mc*sizeof(int));
it_acc_vol = malloc(num_red_mc*sizeof(int));
it_acc_mc = malloc(num_red_mc*sizeof(int));
save_en = dbl_aloc2d(no_save_en,10);
save_each_en = dbl_aloc2d(num_red_mc,4);


order(names,ord);



if(selection == 0){		// Select individual bsc1 stiffness matrices
readsq_ind(&itot,ist,seq_file);

readss_ind(stif,seq1,geom,ord,seq_file,stif_file);

assign_ind(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 1){		// Select average bsc0 dimer stiffness matrices
readsq_di(&itot,ist,seq_file);

readss_di(stif,seq1,geom,ord);

assign_di(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 2){		// Select average bsc0 tetramer stiffness matrices
readsq_tetra(&itot,ist,seq_file);

readss_tetra(stif,seq1,geom,ord);

assign_tetra(itot,ist,seq1,stif,stiff,geom,xconf0);
}


if(selection == 3){		// Select average bsc0 tetramer stiffness matrices with bsc0 dimer at end
readsq_tetra_di(&itot,ist,seq_file);

readss_tetra_di(stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,ord);

assign_tetra_di(itot,ist,stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,stiff,xconf0);
}



// Get stiffness tetramer values
readss_tetra(stif_tetra_gen,seq_tetra_gen,geom_tetra_gen,ord);
// Get coordinates for nucleosome
read_nucl(len_bp_nucl,xconf_nucl,ist_nucl,&len_seq_nucl);

readini(itot,xconfi);
// Set working coordinates and energies equal to the starting ones
transfer_hel_coord(0,itot,xconf,xconfi);
// Prepare MC sequence xconf0_mc, stiff_mc
prepare_mc(xconf0,xconf0_mc,ist_mc,ist,ist_nucl,ins_nucl,no_nucl,len_seq_nucl,geom_tetra_gen,seq_tetra_gen,itot,itot_mc,stiff,stiff_mc,stif_tetra_gen);
// Initialize starting linker configuration
readini(itot_mc,xconfi_mc);
// Set working coordinates and energies equal to the starting ones
transfer_hel_coord(0,itot_mc,xconf_mc,xconfi_mc);


num_bp = itot+1;

// Compute energy
	energy(itot_mc,xconf_mc,xconf0_mc,dx,stiff_mc,ener);


	double ener0=0.0;
	double enertot;

	for(i=0; i<itot_mc; i++){ener0 = ener0 + ener[i];}

	enertot = ener0;

	buf_enertot = enertot;	
	printf("enertot %lf\n",enertot);


	
	printf("num_red_mc %d\n",num_red_mc);

	/* --- Introduce floating particle --- */
	double br_mot[3];	// Initial values floating particle
	br_mot[0] = 0.0;
	br_mot[1] = 0.0;
	br_mot[2] = 0.0;
	
	double rx,ry,rz;	

/* ------------------------------- Looping starts ------------------------------- */

ind_ener = malloc(it_snapshots*sizeof(double));

for(phi=0;phi<it_prog;phi++){			// do optimization it_prog times to check if MC converges


// Initialize variables
//xconf_all, pos_old, orient_old, nucl_pos, R_T_mst,R_T_i,d_nucl,d_nucl_bp,d_nucl_fl,center_bp,d_bp,d_bp_fl
init_dbl2d_arr(xconf_all,6,len_all);
init_dbl2d_arr(pos_old,len_all,3);
init_dbl3d_arr(orient_old,len_all,3,3);
init_dbl2d_arr(nucl_pos,no_nucl,3);
init_dbl3d_arr(R_T_mst,len_all,3,3);
init_dbl3d_arr(R_T_i,len_all,3,3);
init_dbl2d_arr(pos_old,len_all,3);
init_dbl2d_arr(d_nucl,no_nucl,no_nucl);
init_dbl2d_arr(d_nucl_bp,no_nucl,red_bp_all);
init_dbl2d_arr(d_nucl_fl,no_nucl,bp_all);
init_dbl2d_arr(d_nucl_ph,no_nucl,2*red_ph);
init_dbl2d_arr(d_bp,red_bp_all,red_bp_all);
init_dbl2d_arr(d_bp_fl,bp_all,bp_all);
init_dbl2d_arr(center_bp,bp_all,6);

if(incr == 1){acc = (phi+1)*7;	
	   acc_ph = 2*acc;	
	   for(i=0;i<no_nucl;i++){q_nucl[i] = 52/acc;}
	      }


	// Compute energy
	energy(itot_mc,xconf_mc,xconf0_mc,dx,stiff_mc,ener);


	double ener0=0.0;
	double enertot;

	for(i=0; i<itot_mc; i++){ener0 = ener0 + ener[i];}

	enertot = ener0;

	buf_enertot = enertot;	


	/* ---- Generating starting structure ---- */

	if(rnd_start == 3)
	{
	// Set working coordinates and energies equal to the starting ones
	transfer_hel_coord(0,itot_mc,xconf_mc,xconfi_mc);
	nucl_ind = no_nucl-1;
	cmb_link_nucl_seg(xconf_mc,xconf_nucl,xconf_all,ist_mc,ist_nucl,ist_all,ins_nucl,no_nucl,len_bp_nucl,len_seq_nucl,itot_mc, len_all, nucl_ind, &len_seg, ind_linker); // To know right segmentation of ind_linker
	
	int buf_it_mcsteps = 1000; // normally 500000
	int buf_metr_on = 1;	// Metropolis on
	int start,stop,buf_len_seg,bufbuf_len_seg,bufbufbuf_len_seg,count_olap, lim_backset;
	nucl_ind = 0;
	int red_bp_all_seg;
	count_olap = 0;

	lim_backset = 1000;
	buf_len_seg = 0;
	bufbuf_len_seg = 0;
	bufbufbuf_len_seg = 0;


	len_seg = 0;	// For starting the loop

	//for(i=0;i<2*(no_nucl+1);i++){printf("%d ",ind_linker[i]);}
	
	k=0;
	for(nucl_ind = 0;nucl_ind < no_nucl;nucl_ind++)	// Loop for contructing structure
	{

	init_dbl2d_arr(d_nucl,no_nucl,no_nucl);
	init_dbl2d_arr(d_nucl_bp,no_nucl,red_bp_all);
	init_dbl2d_arr(d_nucl_fl,no_nucl,bp_all);
	init_dbl2d_arr(d_nucl_ph,no_nucl,2*red_ph);
	init_dbl2d_arr(d_bp,red_bp_all,red_bp_all);
	init_dbl2d_arr(d_bp_fl,bp_all,bp_all);
	init_dbl2d_arr(center_bp,bp_all,6);

	energy(itot_mc,xconf_mc,xconf0_mc,dx,stiff_mc,ener);
	k = 2*nucl_ind;
	start = ind_linker[k]; 
	stop = ind_linker[k+1]; 
	buf_it_mcsteps = 10*(stop-start);
	if(count_olap > lim_backset/2){buf_it_mcsteps = 100*(stop-start);}
	if(count_olap > lim_backset/1.5){buf_it_mcsteps = 10*(stop-start);}
	transfer_hel_coord_gen(start,stop,xconf_mc,start,xconfi_mc);
	monte_carlo_seg(buf_it_mcsteps,start,stop,enertot,ener,xconft_mc,xconf_mc,sc,xconf0_mc,stiff_mc,enex,dd,boltz,temp, &endenergy, buf_metr_on);

	if(nucl_ind == 1){transfer_hel_coord_gen(start,stop,xconf_mc,start,xconfi_mc);}
	
	//bufbufbuf_len_seg = bufbuf_len_seg;
	//bufbuf_len_seg = buf_len_seg;
	//buf_len_seg = len_seg;
	cmb_link_nucl_seg(xconf_mc,xconf_nucl,xconf_all,ist_mc,ist_nucl,ist_all,ins_nucl,no_nucl,len_bp_nucl,len_seq_nucl,itot_mc, len_all, nucl_ind, &len_seg, ind_linker);

	
	//for(i=0;i<2*(no_nucl+1);i++){printf("%d ",ind_linker[i]);}
	

	cart_rec_old(len_seg,no_rot_pars,xconf_all,pos_old,orient_old,R_T_mst,R_T_i);	

	nucl_positions(nucl_pos,orient_old,pos_old,center_nucl,ins_nucl,nucl_ind+1);
	//for(i=0;i<10;i++){printf("%f %f %f\n", nucl_pos[i][0], nucl_pos[i][1], nucl_pos[i][2]);}
	nucl_dist(len_seg+1,nucl_ind+1,num_fl,len_seg,acc,center_bp,nucl_pos,ins_nucl,pos_old,d_nucl,d_nucl_bp,d_nucl_fl,xyz_float, cut_d_bp_nucl,&m_nucl_bp);


	
	red_bp_all_seg = (int) (double)(len_seg)/(double)acc;

	//printf("nucl_ind %d start %d stop %d len_seg %d count_olap %d red_bp_all_seg %d\n",nucl_ind,start,stop,len_seg,count_olap,red_bp_all_seg);

	check_olap(len_seg+1,num_fl,len_seg,acc,center_bp,pos_old,d_bp,d_bp_fl,xyz_float,bp_excl_vol,bp_fl_excl_vol, nucl_ind,d_nucl,d_nucl_bp,d_nucl_fl,nucl_nucl_excl_vol, nucl_bp_excl_vol,nucl_fl_excl_vol,red_bp_all_seg,&vol_olap_start);

	printf("Overlap pre-starting structure %d\n",vol_olap_start);

	if(vol_olap_start == 0){
				write_bp_pos_in_table_seg(table_bp,len_seg,phi,pos_old,out_folder,nucl_ind);
				write_nucl_pos_in_table_seg(table_nucl,no_nucl,phi,nucl_pos,out_folder,nucl_ind);
				count_olap = 0;
				}	

		if(vol_olap_start == 1){
					//len_seg = buf_len_seg;
					count_olap++;
					if(count_olap > lim_backset){if(nucl_ind > 2){
										nucl_ind = nucl_ind - 2; 
										count_olap = 0; 										printf("nucl_ind %d yes\n", nucl_ind);
											}
								     if(nucl_ind == 2){
										nucl_ind = nucl_ind - 1; 
										count_olap = 0; 										printf("nucl_ind %d yes\n", nucl_ind);
											}
								      if(nucl_ind == 1){
										count_olap = 0;
										printf("nucl_ind %d yes\n", nucl_ind);
											}
									}
					if(nucl_ind > 1){nucl_ind = nucl_ind - 2;}
					else{nucl_ind = 0;}

				write_bp_pos_in_table_seg(table_bp,len_seg,phi,pos_old,out_folder,nucl_ind);
				write_nucl_pos_in_table_seg(table_nucl,nucl_ind+1,phi,nucl_pos,out_folder,nucl_ind);

					}

	}
	

	}	// rnd_start = 3 if

	/* ---- END Generating starting structure ---- */


	/*---------------------------------MC algorithm-----------------------------*/

	/* --- Initialize floating particle --- */

	initial_fl_pos(num_fl,shift_in,br_mot,xyz_float);

	// Charge configuration of floating particles

	for(i=0;i<num_fl;i++){
			q_fl[i] = 1.0;
			      }

	/* ----------- Start algorithm -------------- */

	buf_tot_en = 1000000;	// Set starting energy high so that first structure always gets accepted
	buf_enertot = 1000000;
	err_count = 0;
	acc_vol = 0;
	acc_mc = 0;
	it_en = 0;
	it_each_en = 0;
	start_sim = 1;	// For creating a starting structure
	
	printf("num_red_mc %d\n",num_red_mc);

    for(psi=0;psi<num_red_mc;psi++){
	
	clock_t start = clock(), diff; //Timer

	err_count++;
	//if(err_count%1 == 0){printf("%d\n",err_count); printf("psi %d\n",psi);}	// Warum das?

	if(psi%10000 == 0){printf("psi %d start_sim %d\n",psi,start_sim);}

	if(start_sim == 2){psi = psi-1;}

	//printf("psi %d\n",psi);

	/*for(j=0; j<(itot_mc); j++){
	for(i=0; i<6; i++){
			printf("%lf ",xconf_buf[i][j]);  
		      }	
			printf("\n");
		  }
printf("\n");*/


	transfer_hel_coord(0,itot_mc,xconf_buf,xconf_mc);	// Put helical parameters in buffer
	
	


		/*for(j=0; j<itot_mc; j++){
	for(i=0; i<6; i++){
			printf("%lf ",xconf_mc[i][j]);  
		      }	
			printf("\n");
		  }*/

	

	// Do Monte Carlo
	
	energy(itot_mc,xconf_mc,xconf0_mc,dx,stiff_mc,ener);
	//for(i=0;i<itot_mc;i++){printf("%lf\n",ener[i]);}
	if(psi<2){ener0=0.0;
		for(i=0; i<itot_mc; i++){ener0 = ener0 + ener[i];}
		//printf("ener0 %lf\n",ener0);
		enertot = ener0;
		endenergy = 0.043424*enertot;
		}



	if(psi == 0 && rnd_start ==1)
		    {   int buf_it_mcsteps = 500000; // normally 500000
			int buf_metr_on = 1;
			monte_carlo(buf_it_mcsteps,itot_mc,enertot,ener,xconft_mc,xconf_mc,sc,xconf0_mc,stiff_mc,enex,dd,boltz,temp, &endenergy, buf_metr_on, &no_mc_move);
			no_mc_move = 0;
			enertot = endenergy*23;
			//printf("get_itt %d\n",get_itt);
			it_dna_mc[psi] = get_itt;
		     }
	else{
	// ---- Normal execution
	monte_carlo_ind(it_mcsteps,itot_mc,enertot,ener,xconft_mc,xconf_mc,sc,xconf0_mc,stiff_mc,enex,dd,boltz,temp,&endenergy, metr_on, &get_itt);
	enertot = endenergy*23;
	//printf("get_itt %d\n",get_itt);
	it_dna_mc[psi] = get_itt;
	// ---- End normal execution
	}

	
	

	//endenergy = endenergy*23;	// (factor 0.043424: kcal/mol recalculated in eV) To account for coarse-graining the energy


	// Combine linker and nucleosomal DNA
	cmb_link_nucl(xconf_mc,xconf_nucl,xconf_all,ist_mc,ist_nucl,ist_all,ins_nucl,no_nucl,len_bp_nucl,len_seq_nucl,itot_mc,len_all,get_itt,&real_itt);

	
	if(psi>1){
	end_to_end_dist_fast(len_all,real_itt,no_rot_pars,xconf_all,pos_old,orient_old,R_T_mst,R_T_i);
		
	         }
	else{	// To reconstruct in the beginning
	cart_rec_old(len_all,no_rot_pars,xconf_all,pos_old,orient_old,R_T_mst,R_T_i);
	//end_to_end_dist_fast(len_all,0,no_rot_pars,xconf_all,pos_new,orient_new,R_T_mst,R_T_i); //printf("exec\n");
	transfer_pos_and_orient(0,len_all,orient_store,pos_store,orient_old,pos_old);
	transfer_R(0,len_all,store_R_T_mst,store_R_T_i,R_T_mst,R_T_i);
	}


	nucl_positions(nucl_pos,orient_old,pos_old,center_nucl,ins_nucl,no_nucl);
	nucl_dist(bp_all,no_nucl,num_fl,itot,acc,center_bp,nucl_pos,ins_nucl,pos_old,d_nucl,d_nucl_bp,d_nucl_fl,xyz_float,cut_d_bp_nucl,&m_nucl_bp);
	
	// Check hardcore overlap between nucleosomes, bp and floating particles

	check_olap(bp_all,num_fl,len_all,acc,center_bp,pos_old,d_bp,d_bp_fl,xyz_float,bp_excl_vol,bp_fl_excl_vol, no_nucl,d_nucl,d_nucl_bp,d_nucl_fl,nucl_nucl_excl_vol, nucl_bp_excl_vol,nucl_fl_excl_vol,red_bp_all,&vol_olap_start);


	if(vol_olap_start == 1){if(psi>0){enertot = buf_enertot;}	// rejected (Volume overlaps)
				if(psi==0){start_sim = 2;}
				transfer_hel_coord(0,itot_mc,xconf_mc,xconf_buf);
				transfer_pos_and_orient(real_itt,len_all,orient_old,pos_old,orient_store,pos_store); 
				if(rnd_start != 1){transfer_R(real_itt,real_itt+1,R_T_mst,R_T_i,store_R_T_mst,store_R_T_i);} 
				//printf("Start new\n");
				
				}
     	else{

	//if(psi==0){start_sim = 1;}

	if(psi>0){it_acc_vol[acc_vol] = get_itt;	
		  acc_vol++;}	// Acceptance counter


	/* ---- Separate linker and nucleosomal DNA ---- */

	sep_nucl_dna(linker_pos, linker_orient, nucl_dna_pos, nucl_dna_orient, orient_old, pos_old, ins_nucl, no_nucl, len_all);

	index_nucl_dna(link_nucl_index,ins_nucl,no_nucl,len_all);
	for(i=0;i<2*no_nucl+2;i++){ph_link_nucl_index[i] = 2*link_nucl_index[i];}

	charges_ph(arr_q_ph, q_link, q_nucl_dna, ph_link_nucl_index, ins_nucl, no_nucl, len_all);


	if(psi>1){
			partial_ph_pos(bp_all,len_all,d_ph_midbp,mid_bp,center_bp,dir_bp,orient_old,dir_ph,pos_ph,real_itt);
		  }
	else	 {	ph_pos(bp_all,len_all,d_ph_midbp,mid_bp,center_bp,dir_bp,orient_old,dir_ph,pos_ph);
			transfer_ph_pos(0,2*bp_all,buf_pos_ph,pos_ph);
		  }
	
	if(min_ph == 1){dist_ph_MIN(bp_all,num_fl,acc_ph,pos_ph,d_ph1,d_ph2,d_ph_fl,xyz_float,no_nucl,ins_nucl,d_nucl_ph,nucl_pos,in_q_ph_ph_MIN,in_q_nucl_ph,cut_d_ph_nucl,&m_nucl_ph,red_ph);
			}
	else{		dist_ph(bp_all,num_fl,acc_ph,pos_ph,d_ph1,d_ph2,d_ph_fl,xyz_float,no_nucl,ins_nucl,d_nucl_ph,nucl_pos,in_q_ph_ph,in_q_nucl_ph,cut_d_ph_nucl,&m_nucl_ph,red_ph);
	     }


	calc_lj_nucl(&k5,&k6,&k7,red_bp_all,no_nucl,d_nucl,d_nucl_bp,d_nucl_fl,lj_nucl_cut_low,lj_nucl_cut_up,lj_nucl,lj_nucl_bp,lj_nucl_fl,lj_k_nucl,lj_k_nucl_bp, lj_k_nucl_fl,lj_a_nucl,lj_a_bp_nucl,acc,num_fl,lj_a_nucl_fl,&e_nucl_lj);

	/* ---------- Calculate DH potential (nucleosome vs nucleosome, ph, float) ---------- */

if(min_ph == 1){calc_dh_nucl_MIN(&k8,&k9,&k10,bp_all,no_nucl,d_nucl,d_nucl_ph,d_nucl_fl,dh_nucl_cut_low,dh_nucl_cut_up,dh_nucl,dh_nucl_ph,dh_nucl_fl,q_nucl, arr_q_ph,q_fl,in_q_nucl_ph,acc_ph,num_fl,eps,kappa,f_pi_eps0,red_ph,&e_nucl_dh);
		}

else{	calc_dh_nucl(&k8,&k9,&k10,bp_all,no_nucl,d_nucl,d_nucl_ph,d_nucl_fl,dh_nucl_cut_low,dh_nucl_cut_up,dh_nucl,dh_nucl_ph,dh_nucl_fl,q_nucl, arr_q_ph,q_fl,in_q_nucl_ph,acc_ph,num_fl,eps,kappa,f_pi_eps0,red_ph,&e_nucl_dh);
     }


	/* ---- Calculate Lenard-Jones (bp + float) ---- */
	
	calc_lj(&k1,&k2,bp_all,d_bp,lj_cut_low,lj_cut_up,lj_bp,lj_k,lj_a_bp,acc,num_fl,lj_bp_fl,d_bp_fl,lj_a_fl,red_bp_all,&e_lj_bp_bpfl);

	/* ---------- Calculate DH potential (bp + float) ---------- */

if(min_ph == 1){calc_dh_MIN(&k3,&k4,bp_all,dh_ph,d_ph1,d_ph2,dh_cut_low,dh_cut_up,arr_q_ph,eps,kappa,f_pi_eps0,acc_ph,num_fl,dh_ph_fl,d_ph_fl,q_fl,in_q_ph_ph_MIN,red_ph,&e_dh_ph_phfl);
		}

else{	calc_dh(&k3,&k4,bp_all,dh_ph,d_ph1,d_ph2,dh_cut_low,dh_cut_up,arr_q_ph,eps,kappa,f_pi_eps0,acc_ph,num_fl,dh_ph_fl,d_ph_fl,q_fl,in_q_ph_ph,red_ph,&e_dh_ph_phfl);
     }

	
	double buf_tot_dh_lj;
	buf_tot_dh_lj = tot_dh_lj;

	tot_dh = 0;
	tot_lj = 0;
	tot_dh_lj = 0;

	tot_dh_lj = e_lj_bp_bpfl + e_dh_ph_phfl;
	
	tot_en = tot_dh_lj + endenergy/acc;			// Energy DH+LJ DNA + elastic DNA (in eV)

	tot_en = tot_en + e_nucl_dh + e_nucl_lj;	// Add interactions of nucleosome

	//tot_en = tot_en*23;				// Recalculate total energy in kcal/mol (from eV)

	metropolis(boltz/23,temp,buf_tot_en,tot_en,&iflag);	// Metropolis for total energy (in eV)
	//printf("buf_tot_en %f tot_en %f iflag %d\n",buf_tot_en,tot_en,iflag);
	
	//iflag = 1;

	if(do_sim_ann == 1){sim_ann(tmax,psi,&t,factor,boltz,temp,buf_tot_en,tot_en,&iflag);}

				if(psi<=1){start_sim = 2;}
				transfer_hel_coord(0,itot_mc,xconf_mc,xconf_buf);
				enertot = buf_enertot;
				if(psi>1){transfer_pos_and_orient(real_itt,len_all,orient_old,pos_old,orient_store,pos_store);}
				if(psi>1){transfer_R(real_itt,real_itt+1,R_T_mst,R_T_i,store_R_T_mst,store_R_T_i);}
				if(psi>1){transfer_ph_pos(2*real_itt, 2*bp_all,pos_ph,buf_pos_ph);}
	      		     }
		else{buf_tot_en = tot_en;	// accepted
		     buf_enertot = enertot;
			if(psi<=1){start_sim = 1;}
			if(psi>1){transfer_pos_and_orient(real_itt,len_all,orient_store,pos_store,orient_old,pos_old);}	
			if(psi>1){transfer_R(real_itt,real_itt+1,store_R_T_mst,store_R_T_i,R_T_mst,R_T_i);}
			if(psi>1){transfer_ph_pos(2*real_itt, 2*bp_all,buf_pos_ph,pos_ph);}
			buf_real_itt = real_itt;	// for transfering coordinates		
		     }	// accepted


	
	//printf("Time taken for one MC step MC: %d seconds %d milliseconds\n", msec/1000, msec%1000);

	if(iflag == 1){if(psi%1 == 0){

					}
			err_count = 0;	
			if(psi>0){it_acc_mc[acc_mc] = get_itt;
				  acc_mc++;}
			if(acc_mc%int_en_save == 0 && want_to_save_en == 1){	
							save_en[it_en][0] = psi;
							save_en[it_en][1] = acc_mc; 
							save_en[it_en][2] = e_lj_bp_bpfl; 
							save_en[it_en][3] = e_dh_ph_phfl; 
							save_en[it_en][4] = e_nucl_lj; 
							save_en[it_en][5] = e_nucl_dh;
							save_en[it_en][6] = endenergy/acc; 
							save_en[it_en][7] = tot_en;

							write_bp_pos_in_table(table_bp,len_all,phi,pos_old,out_folder,it_en); 
						//write_bp_triad_in_table(table_bp_triad,len_all,phi,orient_old,out_folder,it_en); 
							write_nucl_pos_in_table(table_nucl,no_nucl,phi,nucl_pos,out_folder,it_en);
							write_helpars_in_table(itot_mc,phi,xconf_mc,out_folder,it_en);
							write_bp_nucl_pos_in_table_GRO(table_pos,len_all,phi,pos_old,no_nucl,nucl_pos,out_folder,it_en);
							//for(j=0;j<8;j++){printf("%lf ",save_en[it_en][j]);} printf("\n");
							it_en++;
						   			    }
			
			if(detect_indep == 1){	save_each_en[it_each_en][0] = psi;
						save_each_en[it_each_en][1] = it_each_en;
						save_each_en[it_each_en][2] = tot_en;
						sed_co = calc_sed_coeff(no_nucl,d_nucl);
						/*inv_d = 0;
						for(i=0;i<no_nucl;i++){
						for(j=i+1;j<no_nucl;j++){
									inv_d = inv_d + 1/d_nucl[i][j];
									}
								       }
						sed_co = 11.1*(1 + 2*54.6/no_nucl * inv_d);*/
					 	save_each_en[it_each_en][3] = sed_co;
						rg = calc_rg(no_nucl,nucl_pos);
						save_each_en[it_each_en][4] = rg;
						//printf("tot_en %f sed_co %f\n",save_each_en[it_each_en][2],save_each_en[it_each_en][3]);
					      	it_each_en++;
						//printf("it_each_en %d\n",it_each_en);
					      }
			}

	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;

	}
	//printf("End if statement vol_olap_start\n");
	
    }


	/* ----- Output positions and triads ----- */
	write_bp_pos_in_table(table_bp,len_all,phi,pos_old,out_folder,-1); 
	write_bp_triad_in_table(table_bp_triad,len_all,phi,orient_old,out_folder,-1); 
	write_nucl_pos_in_table(table_nucl,no_nucl,phi,nucl_pos,out_folder,-1); 
	write_helpars_in_table(itot_mc,phi,xconf_mc,out_folder,-1);
	write_bp_nucl_pos_in_table_GRO(table_pos,len_all,phi,pos_old,no_nucl,nucl_pos,out_folder,-1);
	write_acc_mc_moves_in_table(0,it_dna_mc,num_red_mc,it_mcsteps,it_prog,it_mc,metr_on,do_sim_ann,stif_file,seq_file,out_folder); 
	write_acc_mc_moves_in_table(1,it_acc_vol,acc_vol,it_mcsteps,it_prog,it_mc,metr_on,do_sim_ann,stif_file,seq_file,out_folder);
	write_acc_mc_moves_in_table(2,it_acc_mc,acc_mc,it_mcsteps,it_prog,it_mc,metr_on,do_sim_ann,stif_file,seq_file,out_folder);
	
	printf("it_each_en %d\n",it_each_en);
	if(want_to_save_en == 1){write_energy_in_table(it_en,phi,save_en,out_folder);}
	if(detect_indep == 1){write_indep_in_table(it_each_en,phi,save_each_en,out_folder);}

	/* ----- Output pdb ----- */

	if(selection == 0 || selection == 1){
	output_str_di(output,len_all,ist_all,xconf_all,phi,out_folder);
	}
	if(selection == 2 || selection == 3){
	output_str_tetra(output,itot,ist,xconf,phi,out_folder);
			  }

	//CEHS_build(len_all,ist_all,xconf_all,phi,output,out_folder);


/* ---------- Loop the manipulation part DELETE THIS PART BECAUSE WE ONLY WANT TO GET END STRUCTURES (KEEP CALCULATION OF POTENTIALS!!)---------- */

	count_no_excl = 0;
	it_mod = 0;
	alpha = 0;



/* ---- Write values for output tables ---- */

int ind;

ind=ord[0];						// assigns which hel par is printed in table
write_hel_coord_in_table(table_heli_shif, itot_mc, ord, xconf_mc, phi, ind);

ind=ord[1];
write_hel_coord_in_table(table_heli_slid, itot_mc, ord, xconf_mc, phi, ind);

ind=ord[2];
write_hel_coord_in_table(table_heli_rise, itot_mc, ord, xconf_mc, phi, ind);

ind=ord[3];
write_hel_coord_in_table(table_heli_tilt, itot_mc, ord, xconf_mc, phi, ind);

ind=ord[4];
write_hel_coord_in_table(table_heli_roll, itot_mc, ord, xconf_mc, phi, ind);

ind=ord[5];
write_hel_coord_in_table(table_heli_twis, itot_mc, ord, xconf_mc, phi, ind);




for(i=0;i<itot+1;i++){		//Put strings in file = one line in the file
fputs(table_heli_shif[i],file_shif);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_slid[i],file_slid);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_rise[i],file_rise);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_tilt[i],file_tilt);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_roll[i],file_roll);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_twis[i],file_twis);}





} 
/* ------------------- End of for loop to execute the whole program it_prog times ------------------- */


fclose(file_shif);
fclose(file_slid);
fclose(file_rise);
fclose(file_tilt);
fclose(file_roll);
fclose(file_twis);
//printf("Test1\n");
write_pars_in_file(seq_file,ins_nucl,no_nucl,selection,it_mcsteps,it_prog,it_mc,metr_on,tot_en,num_red_mc,acc_vol,acc_mc,do_sim_ann,tmax,factor,lj_cut_up,dh_cut_up,rnd_start,acc,acc_ph,tot_lj*23,tot_dh*23,out_folder);
//printf("Test2\n");

sum_energy=0.0;
for(i=0;i<it_snapshots;i++){
sum_energy = sum_energy + ind_ener[i];
//printf("%f \n", ind_ener[i]);
}
free (ind_ener);


ind_energy = sum_energy/(double) it_snapshots;


printf("Acceptance rate for Volume overlap: %lf\n", (double) acc_vol/ (double) num_red_mc);
printf("Acceptance rate for overall MC: %lf\n", (double) acc_mc/ (double) num_red_mc);
printf("Time taken for one MC step MC: %d seconds %d milliseconds\n", msec/1000, msec%1000);
//printf("%d steps needed to get %d structures that are %d times energy-accepted by Metropolis = %lf efficiency\n", it_mod, count_no_excl, alpha ,(double) alpha/ (double) it_mod);
printf("Time taken for one of the %d variations with %d MC steps: %d seconds %d milliseconds\n", it_snapshots, it_mcsteps, msec2/1000, msec2%1000);
printf("Time taken for reconstruction: %d seconds %d milliseconds\n", msec1/1000, msec1%1000);
printf("Mean energy after %d stabilization steps and modifying %d times by %d MC steps: %f eV\n", it_mc, it_snapshots, it_mcsteps, ind_energy);



diff_all = clock() - start_all;

int msec_all = diff_all * 1000 / CLOCKS_PER_SEC;
printf("Time taken for whole algorithm: %d seconds %d milliseconds\n", msec_all/1000, msec_all%1000);


/* Free memory of arrays */

free(seq1);
free(ist);
free_dbl3d(stif,alc_space);
free_dbl3d(stiff,len_seq);
free(geom);
free(xconfi);
free(xconf0);
free(xconf);
free(xconft);
free(ener);
free(dx);
free(output);

free(seq1_tetra);
free(seq1_di);
free_dbl3d(stif_tetra,256);
free_dbl3d(stif_di,16);
free(geom_tetra);
free(geom_di);
//printf("Test11\n");
free(table_heli_shif);
free(table_heli_slid);
free(table_heli_rise);
free(table_heli_tilt);
free(table_heli_roll);
free(table_heli_twis);
//printf("Test12\n");
free(out_part);
free(xyz_float);
free(q_fl);
free(center_bp);
/*free(pdb_data);
free(tmp_at);
free(baset);
free(strand);
free(base_num);
free(tmp_xyz);*/
free(d_bp);
free(d_bp_fl);
free(lj_bp);
free(lj_bp_fl);
free(mid_bp);
free(dir_bp);
free(dir_ph);
free(pos_ph);
free(buf_pos_ph);
free(d_ph1);
free(d_ph2);
free(d_ph_fl);
free(dh_ph);
free(dh_ph_fl);
//printf("Test13\n");
free(xconf_buf);
free(xyz_float_buf);

free(pdb_tot);
//printf("Test14\n");
free_dbl3d(orient_old,len_all);
free_dbl3d(orient_new,len_all);
free_dbl3d(orient_store,len_all);
free(pos_old);
free(pos_new);
free(pos_store);
free_dbl3d(rot_new,2*len_all);
//free_dbl3d(rot_store,2*len_all);
free_dbl3d(store_R_T_i,len_all);
free_dbl3d(store_R_T_mst,len_all);
free_dbl3d(R_T_i,len_all);
free_dbl3d(R_T_mst,len_all);
free_dbl3d(buf_R_T_i,len_all);
free_dbl3d(buf_R_T_mst,len_all);
//printf("Test15\n");
free(xconf_nucl);
free(ist_nucl);
free(ist_all);
free(xconf_all);
free(seq_tetra_gen);
free_dbl3d(stif_tetra_gen,256);
free(geom_tetra_gen);
//printf("Test16\n");
free(xconf_mc);
free(xconf0_mc);
free_dbl3d(stiff_mc,itot_mc);
free(xconfi_mc);
free(xconft_mc);
free(ist_mc);

free(nucl_pos);
free(center_nucl);
free(tr_mat);
free(buf_vec);
free(d_nucl);
free(d_nucl_bp);
free(d_nucl_fl);
free(d_nucl_ph);

free(table_bp);
free(table_pos);
free(table_bp_triad);
free(table_nucl);

free(it_dna_mc);
free(it_acc_vol);
free(it_acc_mc);

free_dbl3d(linker_orient,len_all);
free(linker_pos);
free_dbl3d(nucl_dna_orient,len_all);
free(nucl_dna_pos);
free(link_nucl_index);
free(ph_link_nucl_index);
free(arr_q_ph);

free(lj_nucl);
free(lj_nucl_bp);
free(lj_nucl_fl);

free_int3d(in_q_ph_ph,2*red_ph);
free(in_q_ph_ph_MIN);
free(in_q_nucl_ph);
free(q_nucl);
free(ind_linker);
free(dh_nucl);
free(dh_nucl_ph);
free(dh_nucl_fl);

free(save_en);
free(save_each_en);

/* End free memory */



}
