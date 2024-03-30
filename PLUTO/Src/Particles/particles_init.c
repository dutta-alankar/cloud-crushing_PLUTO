/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Initialize particle distrbution, set B.C. and injection.
 
 This file contains routines to initialize particles on the grid,
 assign user-defined boundary conditions and inject particles.
 Particles attributes that can be set are: position, velocity, color.
 In case of evolution with spectra, the initial spectral profile is also
 prescribed here.
 
 \authors A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
          D. Mukherjee
          
 \date    July 24, 2018
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_Init(Data *d, Grid *grid)
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{
  int i,j,k, np, dir, nc;
  int np_glob = RuntimeGet()->Nparticles_glob;
  int np_cell = RuntimeGet()->Nparticles_cell;
  static int first_call = 1;
  double xbeg[3], xend[3];
  Particle p;

  if (first_call){
    RandomSeed(time(NULL),0);
    first_call = 0;
  }

/* --------------------------------------------------------------
   1. Global initialization
   -------------------------------------------------------------- */

  if (np_glob > 0){

    for (dir = 0; dir < 3; dir++){
      xbeg[dir] = grid->xbeg_glob[dir];
      xend[dir] = grid->xend_glob[dir];
    }

    for (np = 0; np < np_glob; np++){
      Particles_LoadUniform(np, np_glob, xbeg, xend, p.coord);
      #if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == YES)
      Particles_LP_InitSpectra(&p);  
      for (nc = 0; nc < PARTICLES_LP_NCOLORS; nc++) p.color[nc] = 0.0;
      #endif
      #if PARTICLES == PARTICLES_CR
      p.speed[IDIR] = 0.0;
      p.speed[JDIR] = 0.0;
      p.speed[KDIR] = 0.0;
      p.mass        = 1.e-3*grid->dV[k][j][i];;
      p.color       = 0.0;
      #endif
      Particles_Insert (&p, d, PARTICLES_CREATE, grid);
    }
  }

/* ------------------------------------------------------------------
   2. Cell by cell initialization.
      Note: You may use Particles_LoadRandom() to initialize
            velocity components but not spatial coordinates.
   ------------------------------------------------------------------ */

  if (np_cell > 0){

    DOM_LOOP(k,j,i){
      xbeg[IDIR] = grid->xl[IDIR][i]; xend[IDIR] = grid->xr[IDIR][i];
      xbeg[JDIR] = grid->xl[JDIR][j]; xend[JDIR] = grid->xr[JDIR][j];
      xbeg[KDIR] = grid->xl[KDIR][k]; xend[KDIR] = grid->xr[KDIR][k];

      for (np = 0; np < np_cell; np++){

      /* -- Spatial distribution -- */

        Particles_LoadUniform(np, np_cell, xbeg, xend, p.coord);

        #if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == YES)
        Particles_LP_InitSpectra(&p);
        for (nc = 0; nc < PARTICLES_LP_NCOLORS; nc++) p.color[nc] = 0.0;
        #endif

     /* -- Velocity distribution -- */

        #if PARTICLES == PARTICLES_CR
        p.speed[IDIR] = RandomNumber(-1,1);
        p.speed[JDIR] = RandomNumber(-1,1);
        p.speed[KDIR] = RandomNumber(-1,1);
        p.rho         = 1.e-3/np_cell;
        #endif
        Particles_Insert (&p, d, PARTICLES_CREATE, grid);
      }
    }
  }  
  Particles_SetID(d->PHead);
}

#if PARTICLES_LP_SPECTRA == YES
/* ********************************************************************** */
void Particles_LP_InitSpectra(Particle* pl)
/*!
 *  Initialize spectra for each particle (only for LAGRANGIAN).
 *  Specify here the initial distribution of N(E) with E for each particle
 *
 *  \param [in]      pl      Pointer to the Particle structure.
 * 
 ********************************************************************** */
{
  int i;
  double Emin, Emax, DeltaE, N_0, alpha1, alpha2, Eb, N1, N2;
  double lnEmin, lnEmax, dlnE, scrh;
  Emin = 1.0e-2;
  Emax = 1.0e4;
  lnEmin = log10(Emin);
  lnEmax = log10(Emax);
  dlnE = (lnEmax - lnEmin)/((double)PARTICLES_LP_NEBINS);
    
  pl->nmicro = 0.001; /* The number density of micro-particles */

/* --------------------------------------------------------
   0. Initialize structure members
   -------------------------------------------------------- */

  pl->cmp_ratio = 1.0;
  pl->shkflag = 0;
  pl->Vshk_upst[RHO] = -1.0;
  pl->Vshk_dnst[RHO] = -1.0;
  pl->cr = 0.0;

/* --------------------------------------------------------
   1. Initialize energy grid.
      Note that pl->eng[i] gives the location of the
      *left* interface of an energy bins.
      There're therefore NBINS+1 nodes in the grid.
   -------------------------------------------------------- */

  for (i = 0; i <= PARTICLES_LP_NEBINS; i++){
    scrh       = lnEmin + i*dlnE;
    pl->eng[i] = pow(10.0, scrh);

    /* Relativistic Maxwellian Distribution*/
    /* double mu = 2.0; // mu = m_0c^2/kT_e 
     double k2_mu = BesselKn(2,mu); 
     pl->chi[i] = (pl->nmicro)*(mu/k2_mu)*(pl->eng[i] + 1.0)
                  *sqrt(pl->eng[i]*(pl->eng[i] + 2.0))*exp(-mu*(pl->eng[i] + 1.0)); 
     pl->cmp_ratio = 1.0; 
     pl->shkflag = 0;
     pl->rho_min = -1.0;
     pl->rho_max = -1.0;
     pl->cr = 0.0; */

  }

/* ----------------------------------------------------------------------
   2. Initialize spectral distribution to power law.
      Chi in each bin is assigned as an average over the energy interval.
      chi_i = \int^Eh_El E^(-p) dE /(Eh - El)
   ---------------------------------------------------------------------- */
  alpha1 = 3.0;
  scrh = (pl->nmicro)/(pow(Emax,1.0-alpha1)-pow(Emin,1.0-alpha1));
  for (i = 0; i < PARTICLES_LP_NEBINS; i++){
    pl->chi[i]  = scrh*(  pow(pl->eng[i+1],1.0-alpha1) 
                        - pow(pl->eng[i],1.0-alpha1))/(pl->eng[i+1] - pl->eng[i]);
  }

}
#endif

/* ********************************************************************* */
void Particles_Inject(Data *data, Grid *grid)
/*!
 *  Inject particles as you wish.
 *
 *  \param [in]  data    Pointer to the PLUTO data structure.
 *  \param [in]  grid    Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Particles_UserDefBoundary(Data *d, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int    dir;
  double xbeg[3], xend[3];
  particleNode *curr = d->PHead, *next;
  Particle *p;
  
  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  if (side == X1_BEG){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X1_END){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_BEG){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_END){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_BEG){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_END){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }
}
