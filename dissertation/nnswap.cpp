/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* umbrella = simple example of how an umbrella program
              can invoke LAMMPS as a library on some subset of procs
   Syntax: umbrella P in.lammps
           P = # of procs to run LAMMPS on
               must be <= # of procs the umbrella code itself runs on
           in.lammps = LAMMPS input script
   See README for compilation instructions */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "library.h"        /* this is a LAMMPS include file */
#include "lammps.h"
#include "math.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)

using namespace LAMMPS_NS;

int main(int narg, char **arg)
{
  int nunit_x = 9;
  int nunit_y = 9;
  int nunit_p = 6;

  /* setup MPI and various communicators
     umbrella is all procs in MPI_COMM_WORLD
    comm_lammps only has 1st P procs (could be all or any subset) */

  MPI_Init(&narg,&arg);

  if (narg != 4) {
    printf("Syntax: umbrella P in.lammps\n");
    exit(1);
  }

  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  int nprocs_lammps = atoi(arg[1]);
  if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int lammps;
  if (me < nprocs_lammps) lammps = 1;
  else lammps = MPI_UNDEFINED;
  MPI_Comm comm_lammps;
  MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);

  /* open LAMMPS input script */

  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[2],"r");
    if (fp == NULL) {
      printf("ERROR: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  FILE *fp_out;
  if (me == 0) {
    fp_out = fopen("my_c_driver.out","w");
    if (fp_out == NULL) {
      printf("ERROR: Could not open out1.txt\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }  

  /* run the input script thru LAMMPS one line at a time until end-of-file
     umbrella proc 0 reads a line, Bcasts it to all procs
     (could just send it to proc 0 of comm_lammps and let it Bcast)
     all LAMMPS procs call lammps_command() on the line */

  void *ptr; 
  void *rnd_ptr;
  if (lammps == 1) lammps_open(0,NULL,comm_lammps,&ptr,&rnd_ptr);
  //LAMMPS *lammps_ptr = (LAMMPS *) ptr;

  int n;
  char line[1024];
  while (1) {
    if (me == 0) {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
      if (n == 0) fclose(fp);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
    if (lammps == 1) lammps_command(ptr,line);
    fprintf(fp_out, line);
  }

  /* run 10 more steps
     get coords from LAMMPS
     change coords of 1st atom
     put coords back into LAMMPS
     run a single step with changed coords */
 
  int* neg;
  int* pos;
  double pe, pe1, temp, temp1, vol, vol1, press, press1, boltz, boltz1, rd, mnm;  
  int *types; 
  double *charges;
  double *coord; 
  int counter = 0;
  int n1,n2,i1,i2,i,oldtneg,oldtpos;
  double oldqneg,oldqpos;
  int unitcell_num, p_num, row, row_up, row_down, col, col_left, col_right;
  int* unitcell_neigh;
  int* p_neigh;

  if (lammps == 1) {
    //lammps_command(ptr,"run 10");

    int natoms = lammps_get_natoms(ptr);

    neg = (int *) malloc(natoms/3*sizeof(int));
    pos = (int *) malloc(2*natoms/3*sizeof(int));
    types = (int *) malloc(natoms*sizeof(int));
    charges = (double *) malloc(natoms*sizeof(double));
    coord = (double *) malloc(3*natoms*sizeof(double));
    unitcell_neigh = (int *) malloc(8*sizeof(int));
    p_neigh = (int *) malloc(6*sizeof(int));

    //if (ptr->logfile) {
//	fprintf(ptr->logfile, "PE!\n");
  //  }
    //if (ptr->screen) {
//	fprintf(ptr->screen, "PE!\n");
  //  }
    //fprintf(fp_out, "natoms: %d\n", natoms);

    while (counter++ < atoi(arg[3])) {
	fprintf(fp_out, "\nCycle %d\n", counter);

	/* record current PE */
    	fprintf(fp_out, "Current PE: ");
    	pe = lammps_get_pe(ptr);
    	fprintf(fp_out, "%f\n", pe);

    	temp = lammps_get_temp(ptr);
	vol = lammps_get_vol(ptr);
	press = lammps_get_press(ptr);

	lammps_get_types(ptr,types);
    	//fprintf(fp_out, "lammps_get_types: done.\n");
    	lammps_get_charges(ptr,charges);
    	//fprintf(fp_out, "lammps_get_charges: done.\n");
    	lammps_get_coords(ptr,coord);
    	//fprintf(fp_out, "lammps_get_coord: done.\n");
    
	/* do the survey of types */
    	n1 = n2 = 0;
    	for (i = 0; i < natoms; i++) {
		if (types[i] == 1) {
			if (charges[i] != -1.4) {
				fprintf(fp_out, "Charge/type inconsistency. Quitting!\n");
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			neg[n1++] = i;
		} else {
			if (charges[i] != 0.7) {
				fprintf(fp_out, "Charge/type inconsistency. Quitting!\n");
				MPI_Abort(MPI_COMM_WORLD,1);
			}			
			pos[n2++] = i;
		}
    	}
    
	/* randomly select a pair of opposite charges */
    	//fprintf(fp_out, "There are %d atoms of type 1 and %d of type 2.\n", n1, n2);
    	i1 = lammps_get_randint(n1, rnd_ptr); 
    	unitcell_num = neg[i1] / nunit_p;
	p_num = neg[i1] % nunit_p;
	row = (int) unitcell_num / nunit_x;
	col = unitcell_num % nunit_x;
	row_up = (row + 1) % nunit_y;
	row_down = (row - 1) % nunit_y;
	col_left = (col - 1) % nunit_x;
	col_right = (col + 1) % nunit_y;
	fprintf(fp_out, "stencil:%d %d %d\n%d %d %d\n", row_down, row, row_up, col_left, col, col_right);
	
	switch (p_num) {
		case 0:
			p_neigh[0] = 	(row_down * nunit_x + col_left)*nunit_p + 6;
			p_neigh[1] =	(row_down * nunit_x + col)*nunit_p + 4;
			p_neigh[2] =	(row * nunit_x + col_left)*nunit_p + 3;
			p_neigh[3] =	(row * nunit_x + col)*nunit_p + 2;
			p_neigh[4] =	(row * nunit_x + col_left)*nunit_p + 6;
			p_neigh[5] =	(row * nunit_x + col)*nunit_p + 4;
			break;

		case 1:
			p_neigh[0] = 	(row_down * nunit_x + col)*nunit_p + 4;
			p_neigh[1] =	(row_down * nunit_x + col)*nunit_p + 5;
			p_neigh[2] =	(row * nunit_x + col)*nunit_p + 1;
			p_neigh[3] =	(row * nunit_x + col)*nunit_p + 3;
			p_neigh[4] =	(row * nunit_x + col)*nunit_p + 4;
			p_neigh[5] =	(row * nunit_x + col)*nunit_p + 5;
			break;

		case 2:
			p_neigh[0] = 	(row_down * nunit_x + col)*nunit_p + 5;
			p_neigh[1] =	(row_down * nunit_x + col)*nunit_p + 6;
			p_neigh[2] =	(row * nunit_x + col)*nunit_p + 2;
			p_neigh[3] =	(row * nunit_x + col_right)*nunit_p + 1;
			p_neigh[4] =	(row * nunit_x + col)*nunit_p + 5;
			p_neigh[5] =	(row * nunit_x + col)*nunit_p + 6;
			break;			

		case 3:
			p_neigh[0] = 	(row * nunit_x + col)*nunit_p + 1;
			p_neigh[1] =	(row * nunit_x + col)*nunit_p + 2;
			p_neigh[2] =	(row * nunit_x + col_left)*nunit_p + 6;
			p_neigh[3] =	(row * nunit_x + col)*nunit_p + 5;
			p_neigh[4] =	(row_up * nunit_x + col)*nunit_p + 1;
			p_neigh[5] =	(row_up * nunit_x + col)*nunit_p + 2;
			break;

		case 4:
			p_neigh[0] = 	(row * nunit_x + col)*nunit_p + 2;
			p_neigh[1] =	(row * nunit_x + col)*nunit_p + 3;
			p_neigh[2] =	(row * nunit_x + col)*nunit_p + 4;
			p_neigh[3] =	(row * nunit_x + col)*nunit_p + 6;
			p_neigh[4] =	(row_up * nunit_x + col)*nunit_p + 2;
			p_neigh[5] =	(row_up * nunit_x + col)*nunit_p + 3;

		case 5:
			p_neigh[0] = 	(row * nunit_x + col)*nunit_p + 3;
			p_neigh[1] =	(row * nunit_x + col_right)*nunit_p + 1;
			p_neigh[2] =	(row * nunit_x + col)*nunit_p + 5;
			p_neigh[3] =	(row * nunit_x + col_right)*nunit_p + 4;
			p_neigh[4] =	(row_up * nunit_x + col)*nunit_p + 3;
			p_neigh[5] =	(row_up * nunit_x + col_right)*nunit_p + 1;
			break;
	}

	n1 = n2 = 0;
	for (i=0; i<6; i++) {
		p_neigh[i] -= 1;
		if (types[p_neigh[i]] == 1) {
			if (charges[p_neigh[i]] != -1.4) {
				fprintf(fp_out, "Charge/type inconsistency. Quitting!\n");
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			n1++;
		} else {
			if (charges[p_neigh[i]] != 0.7) {
				fprintf(fp_out, "Charge/type inconsistency. Quitting!\n");
				MPI_Abort(MPI_COMM_WORLD,1);
			}			
			n2++;
		}
	}

	int * pos_nn = (int *) malloc(n2*sizeof(int));

	n2 = 0;
	for (i=0; i<6; i++) {
		if (types[p_neigh[i]] == 2) {
			pos_nn[n2] = p_neigh[i];		
		}			
		n2++;
	}
	
	i2 = lammps_get_randint(n2, rnd_ptr);
	
	/* find a neighbor of opposite type */
	
	fprintf(fp_out, "Randomly selecting one of type 1: %d (%d)\n", i1, neg[i1]);
    	fprintf(fp_out, "Randomly selecting one of oppositely charged neighbors of %d: %d (%d)\n", neg[i1], i2, pos_nn[i2]);

	/*    
    	for (i = 0; i < natoms; i++) {
      		fprintf(fp_out, "%d, %d, %.1f\n", i, types[i], charges[i]);
    	}
    	*/
	    
    	for (i = 0; i < n2; i++) {
      		fprintf(fp_out, "%d, %d, %.1f\n", i, types[pos_nn[i]], charges[pos_nn[i]]);
    	}
    	
	/* switch the charges */
    	fprintf(fp_out, "Switching types and charges.\n", natoms);
	oldtneg = types[neg[i1]];
	oldqneg = charges[neg[i1]];
    	types[neg[i1]] = 2;
	oldtpos = types[pos_nn[i2]];
	oldqpos = charges[pos_nn[i2]];
    	types[pos_nn[i2]] = 1;
    	charges[neg[i1]] = 0.7;
    	charges[pos_nn[i2]] = -1.4;

	/* write back type/charge info to LAMMPS data-structure */
    	//fprintf(fp_out, "lammps_put_types: start.\n", natoms);
    	lammps_put_types(ptr,types);
    	//fprintf(fp_out, "lammps_put_types: done.\n", natoms);

    	//fprintf(fp_out, "lammps_put_charges: start.\n", natoms);
    	lammps_put_charges(ptr,charges);
    	//fprintf(fp_out, "lammps_put_charges: done.\n", natoms);

	/*
    	for (int i = 0; i < natoms; i++) {
      		fprintf(fp_out, "%d, %d, %.1f\n", i, types[i], charges[i]);
    	}
	*/

	/* compute pe and the resulting boltzmann factor */
    	//lammps_command(ptr, "run 2500");
    	lammps_command(ptr, "run 0");

	/* consistency check */
	lammps_get_types(ptr,types);
   	lammps_get_charges(ptr,charges);
	fprintf(fp_out, "Old type for particle %d was %d, now %d\n", neg[i1], oldtneg, types[neg[i1]]);
	fprintf(fp_out, "Old charge for particle %d was %.1f, now %.1f\n", neg[i1], oldqneg, charges[neg[i1]]);
	fprintf(fp_out, "Old type for particle %d was %d, now %d\n", pos_nn[i2], oldtpos, types[pos_nn[i2]]);
	fprintf(fp_out, "Old charge for particle %d was %.1f, now %.1f\n", pos_nn[i2], oldqpos, charges[pos_nn[i2]]);
    	fprintf(fp_out, "New PE: ");
    	pe1 = lammps_get_pe(ptr);
    	fprintf(fp_out, "%f\n", pe1);

    	temp1 = lammps_get_temp(ptr);
	vol1 = lammps_get_vol(ptr);
	press1 = lammps_get_press(ptr);
    	boltz = exp(-natoms*(pe1-pe)/temp);
	boltz1 = exp(-natoms*(pe1-pe)/temp + natoms*press*(vol1-vol) - natoms*log(vol1/vol)*temp);
	rd = lammps_get_rand(rnd_ptr);
	mnm = MIN( boltz, 1.0 );
    	fprintf(fp_out, "Temperature was %g, volume was %g, pressure was %g, Boltzmann factor1 was %g, Boltzmann factor2 was %g, random number generator gave %g\n", temp, vol, press, boltz, boltz1, rd);

	/* accept or reject the move based on the bolzmann factor */
    	if (rd < mnm) {
		/* accept the move */
    		fprintf(fp_out, "Accepted the move. Relax for 10 steps...\n", natoms);	
    		//fprintf(fp_out, "Accepted the move.\n", natoms);	

		/* relax for 100 */
    		lammps_command(ptr,"run 10");

    	} else {
		/* reject the move, switch types and charges back */
    		fprintf(fp_out, "Rejected the move. Switch types and charges back.\n");
    		types[neg[i1]] = 1;
    		types[pos_nn[i2]] = 2;
    		charges[neg[i1]] = -1.4;
    		charges[pos_nn[i2]] = 0.7;
    		lammps_put_types(ptr,types);
    		lammps_put_charges(ptr,charges);

		/* switch back to old coordinates */
		//lammps_put_coords(ptr,coord);

		/* consistency check */	
    		lammps_get_types(ptr,types);
    		lammps_get_charges(ptr,charges);
		fprintf(fp_out, "Type of particle %d is now %d\n", neg[i1], types[neg[i1]]);
		fprintf(fp_out, "Charge on particle %d is now %.1f\n", neg[i1], charges[neg[i1]]);
		fprintf(fp_out, "Type of particle %d is now %d\n", pos_nn[i2], types[pos_nn[i2]]);
		fprintf(fp_out, "Charge on particle %d is now %.1f\n", pos_nn[i2], charges[pos_nn[i2]]);

		/* relax for 1000 */
		fprintf(fp_out, "Relax for 10 steps...\n");
                lammps_command(ptr,"run 10");
	}

	free(pos_nn);

    vol = lammps_get_vol(ptr);
    press = lammps_get_press(ptr);
    fprintf(fp_out, "Volume is now %g, pressure is now %g\n", vol, press);

    /* end cycling */
    }

    free(types);
    free(charges);
    free(unitcell_neigh);
    free(p_neigh);
  }

  //printf("PE\n");
  //if (lammps->logfile) {
	//fprintf(lammps->logfile, "PE!\n");
  //}

  if (lammps == 1) lammps_close(ptr,rnd_ptr);

  /* close down MPI */

  MPI_Finalize();
  fclose(fp_out);
}
