/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */
#include <stdio.h>      /* printf */
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "compute_contact_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"



#include <iostream>
#include <fstream>
#include<sstream>
#include<string>
using namespace std;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeContactAtom::ComputeContactAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  
  if (narg < iarg) error->all(FLERR,"Illegal compute contact/atom command");

  skin = 0.;

  if(narg > iarg)
  {
      if (narg < iarg+2)
          error->all(FLERR,"Illegal compute contact/atom command");
      if(strcmp("skin",arg[iarg++]))
          error->all(FLERR,"Illegal compute contact/atom command, expecting keyword 'skin'");
      skin = atof(arg[iarg++]);
  }
  
  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 1;

  nmax = 0;
  contact = NULL;

  // error checks

  if (!atom->sphere_flag && !atom->superquadric_flag)
    error->all(FLERR,"Compute contact/atom requires atom style sphere or atom style superquadric!");
}

/* ---------------------------------------------------------------------- */

ComputeContactAtom::~ComputeContactAtom()
{
  memory->destroy(contact);
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute contact/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"contact/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute contact/atom");

  // need an occasional neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
  
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,radsumsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow contact array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(contact);
    nmax = atom->nmax;
    memory->create(contact,nmax,"contact/atom:contact");
    vector_atom = contact;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  // compute number of contacts for each atom in group
  // contact if distance <= sum of radii
  // tally for both I and J

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *tag = atom->tag;
  int ntimestep = update->ntimestep;


  for (int i=0;i<20;i++) rad[i]=100; //it can help sorting


  //radius_sort(rad, nlocal); //calling this function helps obtaining and sorting radius information


  double rad[]={0.0013365, 0.0016839, 0.0019276, 0.0021216, 0.0022854, 0.0024286, 0.0025566, 0.002673, 0.00278, 0.0028794, 0.0029724, 0.0030598, 0.0031426, 0.0032212, 0.0032961, 0.0033678, 0.0034365, 0.0035026, 0.0035663, 0.0036278};

  /* contact_search is for current particle-particle contact binary matrix , initially put to zero */
  for (int i=1;i< nlocal+1;i++)
    for (int j=1; j< 21; j++)
	{
	contact_search[i][j]=0;
	}
 
  /*updating contact_search matrix */
  for (i = 0; i < nall; i++) contact[i] = 0.0; //inbuilt ligght variable

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = radi + radius[j];         
	radsumsq = radsum*radsum;
        if (rsq <= radsumsq) {
          contact[i] += 1.0;
          contact[j] += 1.0;
	  //contact_search[i][counter]=j;
	  //contact_search[j][counter]=i;
        }
      }
    }
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=0; i< nlocal; i++){
  	int counter_search=1;
  	for (int j=0; j< nlocal; j++)
  	{
  	if (i!=j){
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = radius[i] + radius[j];
	radsumsq = radsum*radsum;
          if (rsq <= radsumsq) {
          contact_search[tag[i]][counter_search]=tag[j];
	  counter_search=counter_search+1;
	  }
	  }
        }     
   }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=1;i<=nlocal;i++){
	for (int j=1; j<=20; j++)
	{
		int search_counter=0;
		for (int k=1; k<=20;k++)
		{
		if (contact_search[i][j]==contact_search_old[i][k]) 
		search_counter=1;
		}
		if (search_counter!=1)
		{
		  int i1=tag_to_serial(i), j1=tag_to_serial(contact_search[i][j]), i2, j2;
		  double r1, r2;
		  r1=radius[i1];
		  r2=radius[j1];
		  i2=radiusToNumber(r1, rad);
		  j2=radiusToNumber(r2, rad);
		  N_col_final[i2][j2]+=1;	
		}
		}
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=0; i< 20; i++) {
  for (int j=0; j< 20; j++)
	{
	printf("%d ", N_col_final[i][j]);
	}
  printf("\n");
  }

	//keeping a copy of contact_search as contact_search_old
  for (int i=1; i<=nlocal;i++)
    for (int j=1; j<= 20; j++)
    {
    contact_search_old[i][j]=contact_search[i][j];
    }
    

  int data1=0, data2=1; 
  if (ntimestep==200000){
  	stringstream ss;  
  	ss<< data1;  
   	string fileName;
  	ss>>fileName;
	
    fileName += ".txt"; // important to create .txt file.
    ofstream createFile;
    createFile.open(fileName.c_str(), ios::app);
    for (int i = 0; i < 20; i++)
	{     
	for (int j=0; j< 20;j++)
		{  
			createFile <<N_col_final[i][j]<<" ";
		} 
	createFile <<"\n";
	}
   	}


  if (ntimestep==200000){
  	stringstream ss;  
  	ss<< data2;  
   	string fileName;
  	ss>>fileName;
	
    fileName += ".txt"; // important to create .txt file.
    ofstream createFile;
    createFile.open(fileName.c_str(), ios::app);
    double particle_volume;
    particle_volume = 4.0*3.14/3.0 * radius[0]*radius[0]*radius[0];
    createFile <<nlocal;
   	}

  
  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm_compute(this);
}

/* ---------------------------------------------------------------------- */

int ComputeContactAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    buf[m++] = contact[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    contact[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeContactAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------*/

int ComputeContactAtom::tag_to_serial(int n)
{
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  for (int i=0; i< nlocal; i++)
  {
	if (n==tag[i])
	return i;
  }

}




/* ----------------------------------------------------------------------*/

int ComputeContactAtom::radiusToNumber(double r1, double *rad)
{
  for (int i=0;i<20;i++)
	if (r1==rad[i])
	return i;

}

/* ----------------------------------------------------------------------*/

void ComputeContactAtom::radius_sort(double *rad, int n)
{
  double *radius = atom->radius;
  int counter1=1;
  rad[0]=radius[0];
	//storing different radius values in rad
  for (int i=1; i< n; i++)
  {
	int counter2=0;
	for (int j=0;j<20;j++)
	{
	  if(radius[i]==rad[j]) counter2++;
	}
	if (counter2==0)
	{
	  rad[counter1]=radius[i];
	  counter1++;
	}
  }	
  double temp;
	//sorting
  for(int i=0;i< 20;i++)
	{
		for(int j=i+1;j< 20;j++)
		{
			if(rad[i]>rad[j])
			{
				temp 	=rad[i];
				rad[i]	=rad[j];
				rad[j]	=temp;
			}
		}
	}
}





