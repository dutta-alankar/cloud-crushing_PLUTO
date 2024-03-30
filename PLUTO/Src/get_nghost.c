/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Return the number of ghost zones.

  Return the number of ghost zones employed by the selected numerical
  algorithm. 
  The minimum number for a 2nd-order algorithm is 2.
  Higher-order interpolation scheme may require more zones.
  
  \authors A. Mignone (mignone@to.infn.it)
  \date    Aug 25, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int GetNghost (void)
/*! 
 * Compute the number of ghost zones, depending on the selected
 * scheme.
 *
 *********************************************************************** */
{
  int nghost = 2;

/* --------------------------------------------------------
   Choose stencil based on reconstruction.
   -------------------------------------------------------- */

  #if    (RECONSTRUCTION == MP5)    || (RECONSTRUCTION == PARABOLIC) \
      || (RECONSTRUCTION == MP5_FD) || (RECONSTRUCTION == WENOZ_FD)  \
      || (RECONSTRUCTION == LINEAR && LIMITER == FOURTH_ORDER_LIM) \
      || (RING_AVERAGE_REC > 2)
  nghost = 3;
  #endif

  #if SHOCK_FLATTENING == ONED

   nghost = MAX(4, nghost);

  #elif SHOCK_FLATTENING == MULTID

/* --------------------------------------------------------
    The MULTID shock flattening only need 2 ghost zones.
    However for axisymmetric simulations with CTU 
    3 zones will ensure that flag[][][] will remain 
    symmetric around the axis.
   -------------------------------------------------------- */

   nghost = MAX(3, nghost);
  #endif

/* --------------------------------------------------------
    The following should operate on the static grid 
    version of the code. Add an extra row of boundary
    zones if CTU+CT is selected.
    At least 3 ghost zones.
   -------------------------------------------------------- */

  #ifdef CTU
   #ifdef STAGGERED_MHD
    nghost++;
   #endif
  #endif

/* --------------------------------------------------------
    FARGO PPM needs at least 3 ghost zones 
   -------------------------------------------------------- */

  #ifdef FARGO
   #if FARGO_ORDER == 3
    nghost = MAX(3,nghost);
   #endif
  #endif

/* --------------------------------------------------------
    Particles are requires to propagate for no more than 2
    zones. However, with option #2 of the CR predictor,
    a single particle is advanced for dt*v (instead of
    dt*v/2) thus potentially ending up in the very last
    ghost zones. In this case, the FORCE cannot be properly
    computed, unless an extra zone is added.
   -------------------------------------------------------- */

  #if (PARTICLES == PARTICLES_CR)
  #if PARTICLES_CR_PREDICTOR == 2
  nghost = MAX(3, nghost);
  #endif
  #endif
  
  return nghost;
}

