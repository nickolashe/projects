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

  char buffer[50];
  sprintf(buffer, "mcswap_%s.out", arg[2]);
  FILE *fp_out;
  if (me == 0) {
    fp_out = fopen(buffer,"w");
    if (fp_out == NULL) {
      printf("ERROR: Could not open %s for writing\n", buffer);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }  

  /* run the input script thru LAMMPS one line at a time until end-of-file
     umbrella proc 0 reads a line, Bcasts it to all procs
     (could just send it to proc 0 of comm_lammps and let it Bcast)
     all LAMMPS procs call lammps_command() on the line */

  void *ptr; 
  void *rnd_ptr;
  sprintf(buffer, "lammps_open_%s.out", arg[2]);
  if (lammps == 1) lammps_open(0,NULL,comm_lammps,&ptr,&rnd_ptr, buffer);
  //LAMMPS *lammps_ptr = (LAMMPS *) ptr;

  int n;
  char line[1024];
  while (1) {
    if (me == 0) {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else { n = strlen(line) + 1; fprintf(fp_out, line); }
      if (n == 0) fclose(fp);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
    if (lammps == 1) lammps_command(ptr,line);
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

  if (lammps == 1) {
    //lammps_command(ptr,"run 10");

    int natoms = lammps_get_natoms(ptr);

    neg = (int *) malloc(natoms/3*sizeof(int));
    pos = (int *) malloc(2*natoms/3*sizeof(int));
    types = (int *) malloc(natoms*sizeof(int));
    charges = (double *) malloc(natoms*sizeof(double));
    //coord = (double *) malloc(3*natoms*sizeof(double));

    //if (ptr->logfile) {
//	fprintf(ptr->logfile, "PE!\n");
  //  }
    //if (ptr->screen) {
//	fprintf(ptr->screen, "PE!\n");
  //  }
    //fprintf(fp_out, "natoms: %d\n", natoms);
   /*
   FILE *lib_out;
   sprintf(buffer, "lammps_neighb_%s.out", arg[2]);
    lib_out = fopen(buffer,"w");
    if (lib_out == NULL) {
      printf("ERROR: Could not open %s for writing\n", buffer);
      MPI_Abort(MPI_COMM_WORLD,1);
     }
   fprintf(fp_out, "Starting\n");
   fclose(lib_out);
   */

    while (counter++ < atoi(arg[3])) {
	if (me == 0) {
		fprintf(fp_out, "\nCycle %d\n", counter);

		fprintf(fp_out, "Getting current temp.\n");
        	temp = lammps_get_temp(ptr);
		fprintf(fp_out, "Getting current vol.\n");
        	vol = lammps_get_vol(ptr);
		fprintf(fp_out, "Getting current press.\n");
        	press = lammps_get_press(ptr);
        	fprintf(fp_out, "Temperature was %g, volume was %g, pressure was %g\n", temp, vol, press);
	
		/* record current PE */
		fprintf(fp_out, "Current PE: ");
		pe = lammps_get_pe(ptr);
		fprintf(fp_out, "%f\n", pe);
	}


        /* get charges and types */
	lammps_get_types(ptr,types);
    	//fprintf(fp_out, "lammps_get_types: done.\n");
    	lammps_get_charges(ptr,charges);
    	//fprintf(fp_out, "lammps_get_charges: done.\n");
    	//lammps_get_coords(ptr,coord);
    	//fprintf(fp_out, "lammps_get_coord: done.\n");
    
	/* do the survey of types */
    	n1 = n2 = 0;
    	for (i = 0; i < natoms; i++) {
		if (types[i] == 1) {
			if (charges[i] != -1.4) {
				if (me == 0) {
					fprintf(fp_out, "Charge/type inconsistency. Quitting!\n");
				}
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			neg[n1++] = i;
		} else {
			if (charges[i] != 0.7) {
				if (me == 0) {
					fprintf(fp_out, "Charge/type inconsistency. Quitting!\n");
				}
				MPI_Abort(MPI_COMM_WORLD,1);
			}			
			pos[n2++] = i;
		}
    	}
    
	if (me == 0) {
        	fprintf(fp_out, "There are %d atoms of type 1 and %d of type 2.\n", n1, n2);

		/* randomly select a pair of opposite charges */
		i1 = lammps_get_randint(n1, rnd_ptr); 
		i2 = lammps_get_randint(n2, rnd_ptr); 

		fprintf(fp_out, "Randomly selecting one of type 1: %d (%d)\n", i1, neg[i1]);
		fprintf(fp_out, "Randomly selecting one of type 2: %d (%d)\n", i2, pos[i2]);
    	}

	MPI_Bcast(&i1,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&i2,1,MPI_INT,0,MPI_COMM_WORLD);

	/*    
    	for (i = 0; i < natoms; i++) {
      		fprintf(fp_out, "%d, %d, %.1f\n", i, types[i], charges[i]);
    	}
    	*/

        /* switch the types and charges on i1 and i2 */
        //lammps_get_types(ptr,types);
        //fprintf(fp_out, "lammps_get_types: done.\n");
        //lammps_get_charges(ptr,charges);
        //fprintf(fp_out, "lammps_get_charges: done.\n");
        //lammps_get_coords(ptr,coord);
        //fprintf(fp_out, "lammps_get_coord: done.\n");

	if (me == 0) {
		/* switch the charges */
    		fprintf(fp_out, "Switching types and charges.\n", natoms);
	}

	oldtneg = types[neg[i1]];
	oldqneg = charges[neg[i1]];
	oldtpos = types[pos[i2]];
	oldqpos = charges[pos[i2]];
    	types[neg[i1]] = oldtpos;
    	types[pos[i2]] = oldtneg;
    	charges[neg[i1]] = oldqpos;
    	charges[pos[i2]] = oldqneg;

	/* write back type/charge info to LAMMPS data-structure */
    	//fprintf(fp_out, "lammps_put_types: start.\n", natoms);
        sprintf(buffer, "%s", arg[2]);
    	lammps_put_types(ptr,types,buffer);
    	//fprintf(fp_out, "lammps_put_types: done.\n", natoms);

    	//fprintf(fp_out, "lammps_put_charges: start.\n", natoms);
        sprintf(buffer, "%s", arg[2]);
    	lammps_put_charges(ptr,charges,buffer);
    	//fprintf(fp_out, "lammps_put_charges: done.\n", natoms);

	/*
    	for (int i = 0; i < natoms; i++) {
      		fprintf(fp_out, "%d, %d, %.1f\n", i, types[i], charges[i]);
    	}
	*/

	/* compute pe and the resulting boltzmann factor */
	if (me == 0) {
		fprintf(fp_out, "Relaxing for 0...\n", natoms);
	}
    	lammps_command(ptr, "run 0");

	/* consistency check to check if charges and types have been swapped*/
	lammps_get_types(ptr,types);
   	lammps_get_charges(ptr,charges);
	if (me == 0) {
		fprintf(fp_out, "Old type for particle %d was %d, now %d\n", neg[i1], oldtneg, types[neg[i1]]);
		fprintf(fp_out, "Old charge for particle %d was %.1f, now %.1f\n", neg[i1], oldqneg, charges[neg[i1]]);
		fprintf(fp_out, "Old type for particle %d was %d, now %d\n", pos[i2], oldtpos, types[pos[i2]]);
		fprintf(fp_out, "Old charge for particle %d was %.1f, now %.1f\n", pos[i2], oldqpos, charges[pos[i2]]);
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
		fprintf(fp_out, "Energy difference: %g\n", natoms*(pe1-pe));
		fprintf(fp_out, "Temp. was %g, volume was %g, press. was %g, Boltzmann factor was %g, Boltzmann factor2 was %g, rand. num. gen. gave %g\n", temp1, vol1, press1, boltz, boltz1, rd);
	}

	MPI_Bcast(&rd,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&mnm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	/* accept or reject the move based on the bolzmann factor */
    	if (rd < mnm) {
		if (me == 0) {
			/* accept the move */
    			fprintf(fp_out, "Accepted the move. Relax for 10 steps...\n", natoms);	
		}

		/* relax for 100 */
    		lammps_command(ptr,"run 100");

    	} else {
		if (me == 0) {
			/* reject the move, switch types and charges back */
    			fprintf(fp_out, "Rejected the move. Switch types and charges back.\n");
		}

    		types[neg[i1]] = oldtneg;
    		types[pos[i2]] = oldtpos;
    		charges[neg[i1]] = oldqneg;
    		charges[pos[i2]] = oldqpos;
                sprintf(buffer, "%s", arg[2]);
    		lammps_put_types(ptr,types,buffer);
                sprintf(buffer, "%s", arg[2]);
    		lammps_put_charges(ptr,charges,buffer);

		/* switch back to old coordinates */
		//lammps_put_coords(ptr,coord);

		/* consistency check */	
    		lammps_get_types(ptr,types);
    		lammps_get_charges(ptr,charges);

		if (me == 0) {
			fprintf(fp_out, "Type of particle %d is now %d\n", neg[i1], types[neg[i1]]);
			fprintf(fp_out, "Charge on particle %d is now %.1f\n", neg[i1], charges[neg[i1]]);
			fprintf(fp_out, "Type of particle %d is now %d\n", pos[i2], types[pos[i2]]);
			fprintf(fp_out, "Charge on particle %d is now %.1f\n", pos[i2], charges[pos[i2]]);
		}

		/* relax for 10 */
		fprintf(fp_out, "Relax for 10 steps...\n");
                lammps_command(ptr,"run 0");
	}

    //vol = lammps_get_vol(ptr);
    //press = lammps_get_press(ptr);
    //fprintf(fp_out, "Volume is now %g, pressure is now %g\n", vol, press);

    /* end cycling */
    }

    free(types);
    free(charges); 
    free(pos);
    free(neg);
    //free(coord);
  }

  //printf("PE\n");
  //if (lammps->logfile) {
	//fprintf(lammps->logfile, "PE!\n");
  //}

  if (lammps == 1) lammps_close(ptr,rnd_ptr);

  /* close down MPI */

  MPI_Finalize();
  
  if (me == 0) {
  	fclose(fp_out);
  }
}
