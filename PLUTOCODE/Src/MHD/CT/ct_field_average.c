/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Average staggered magnetic field to zone center.

  \author A. Mignone (mignone@to.infn.it)
  \date   Sep 04, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CT_AverageStaggeredFields (double ****Vs, double ****UU, RBox *box,
                               Grid *grid)
/*!
 * Average staggered magnetic field components to form a zone-centered 
 * field, e.g., \f$ \av{B_i} = (B_{i+1/2} + B_{i-1/2})/2\f$.
 * In cylindrical coordinates the volume-averaged radial component 
 * is obtained as
 * \f[
 *   \av{B_R} = \frac{1}{\Delta V} \int B(R) R dR 
 *    \qquad\mathrm{where}\qquad
 *     B(R) =   B_{R,i+\HALF} \frac{R - R_{i-\HALF}}{\Delta R} 
 *            + B_{R,i-\HALF} \frac{R_{i+\HALF} - R}{\Delta R}
 * \f] 
 * this yields
 * \f[
 *   \av{B_R} =  \frac{R_{i+\HALF} + 2R_{i-\HALF}}{3(R_{i+\HALF}+R_{i-\HALF})}
 *                B_{R,i-\HALF}
 *             + \frac{2R_{i+\HALF} + R_{i-\HALF}}{3(R_{i+\HALF}+R_{i-\HALF})}
 *                B_{R,i-\HALF}
 * \f]
 * Simlar expressions hold in spherical coordinates. 
 * Check with the following MAPLE script:
 * \code
   restart;
   J     := sin(xi); # J = xi, xi^2, sin(xi)
   Br    := (xi - a)/Delta*B[R] + (b - xi)/Delta*B[L];
   Bav   := int (Br*J, xi=a..b) / int (J, xi=a..b):
   Delta := b-a; 
   Bav   := simplify(Bav);
   cp := coeff(Bav, B[R]);
   cm := coeff(Bav, B[L]);
   
   # Check symmetry propetries (sign must not change)
    
   eval([cm,cp],{a=0,b=1}); 
   eval([cm,cp],{a=-1,b=0});
 * \endcode
 * 
 * The averaging is done inside an arbitrary rectangular box.
 * The box may include boundary cells only during the predictor 
 * step of the CT-CTU algorithm, whereas is useless for RK time stepping.
 *
 * When the CT_EN_CORRECTION flag is enabled, we also redefine the
 * zone total energy using the newly formed cell-centered field, i.e.,
 *
 *  \f[ 
 *   E_i \to E_i - \frac{\mathbf{B}_i^2}{2} + \frac{<\mathbf{B}_i^2>}{2}
 *  \f]
 *
 * \param [in]   Vs    array of staggered fields
 * \param [out]  UU    array of conservative variables
 * \param [in]   box   pointer to RBox structure 
 * \param [in]   grid  pointer to Grid structure
 ************************************************************************ */
{
  int i, j, k;
  double b2_old, b2_new;
  double Bx1_ave, Bx2_ave, Bx3_ave;
  double cp, cm, rp, rm;
  double *x1r = grid->xr[IDIR], *x1l = grid->xl[IDIR];
  #if GEOMETRY == SPHERICAL
  double dp, dm, thp, thm;
  double *x2r = grid->xr[JDIR], *x2l = grid->xl[JDIR];
  #endif
  DIM_EXPAND(double ***Bx1s = Vs[BX1s];  ,
             double ***Bx2s = Vs[BX2s];  ,
             double ***Bx3s = Vs[BX3s]; )
  #if PHYSICS == ResRMHD
  double Ex1_ave, Ex2_ave, Ex3_ave;
  DIM_EXPAND(double ***Ex1s = Vs[EX1s];   , 
             double ***Ex2s = Vs[EX2s];   ,
             double ***Ex3s = Vs[EX3s];)
  #endif  
  
/* ---------------------------------------------------------
    Loop over all zones in the given box.
   --------------------------------------------------------- */

  BOX_LOOP(box,k,j,i){
    #if GEOMETRY == CARTESIAN 
    DIM_EXPAND(Bx1_ave = 0.5*(Bx1s[k][j][i] + Bx1s[k][j][i-1]);  ,
               Bx2_ave = 0.5*(Bx2s[k][j][i] + Bx2s[k][j-1][i]);  ,
               Bx3_ave = 0.5*(Bx3s[k][j][i] + Bx3s[k-1][j][i]); )

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(Ex1_ave = 0.5*(Ex1s[k][j][i] + Ex1s[k][j][i-1]);  ,
               Ex2_ave = 0.5*(Ex2s[k][j][i] + Ex2s[k][j-1][i]);  ,
               Ex3_ave = 0.5*(Ex3s[k][j][i] + Ex3s[k-1][j][i]); )
    #endif

    #elif GEOMETRY == CYLINDRICAL

    rp = x1r[i];
    rm = x1l[i];
    cm = (rp + 2.0*rm)/(3.0*(rp + rm));
    cp = 1.0 - cm; 
    Bx1_ave  = cm*Bx1s[k][j][i-1] + cp*Bx1s[k][j][i];
    Bx2_ave  = 0.5*(Bx2s[k][j-1][i] + Bx2s[k][j][i]);

    #elif GEOMETRY == POLAR

    rp = x1r[i];
    rm = x1l[i];
    cm = (rp + 2.0*rm)/(3.0*(rp + rm));
    cp = 1.0 - cm;

    DIM_EXPAND(Bx1_ave = cm*Bx1s[k][j][i-1] + cp*Bx1s[k][j][i];   ,
               Bx2_ave = 0.5*(Bx2s[k][j-1][i] + Bx2s[k][j][i]);   ,
               Bx3_ave = 0.5*(Bx3s[k-1][j][i] + Bx3s[k][j][i]);)

    #elif GEOMETRY == SPHERICAL
    rp = x1r[i];
    rm = x1l[i];

    cm =  (3.0*rm*rm + 2.0*rm*rp + rp*rp)
         /(4.0*(rm*rm + rm*rp + rp*rp));
    cp = 1.0 - cm;
    
    thp = x2r[j];
    thm = x2l[j];
 
    dm  = (thm - thp)*cos(thm) - sin(thm) + sin(thp);
    dm /= (thm - thp)*(cos(thm) - cos(thp));
    dp  = 1.0 - dm;
    DIM_EXPAND(Bx1_ave  = cm*Bx1s[k][j][i-1] + cp*Bx1s[k][j][i];  ,
               Bx2_ave  = dm*Bx2s[k][j-1][i] + dp*Bx2s[k][j][i];  ,
               Bx3_ave  = 0.5*(Bx3s[k-1][j][i] + Bx3s[k][j][i]);)
    #endif
   
   /* ---------------------------------------------
          apply energy correction if necessary 
       --------------------------------------------- */
 
    #if CT_EN_CORRECTION == YES && HAVE_ENERGY
    b2_old = DIM_EXPAND(  UU[k][j][i][BX1]*UU[k][j][i][BX1], 
                        + UU[k][j][i][BX2]*UU[k][j][i][BX2],  
                        + UU[k][j][i][BX3]*UU[k][j][i][BX3]);
    #if PHYSICS == ResRMHD
    b2_old += DIM_EXPAND(  UU[k][j][i][EX1]*UU[k][j][i][EX1], 
                         + UU[k][j][i][EX2]*UU[k][j][i][EX2],  
                         + UU[k][j][i][EX3]*UU[k][j][i][EX3]);
    #endif
    #endif   

    DIM_EXPAND(UU[k][j][i][BX1] = Bx1_ave;  ,
               UU[k][j][i][BX2] = Bx2_ave;  ,
               UU[k][j][i][BX3] = Bx3_ave; )

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(UU[k][j][i][EX1] = Ex1_ave;  ,
               UU[k][j][i][EX2] = Ex2_ave;  ,
               UU[k][j][i][EX3] = Ex3_ave; )
    #endif

    #if CT_EN_CORRECTION == YES && HAVE_ENERGY
    b2_new = DIM_EXPAND(  Bx1_ave*Bx1_ave, 
                        + Bx2_ave*Bx2_ave, 
                        + Bx3_ave*Bx3_ave);
    #if PHYSICS == ResRMHD
    b2_new += DIM_EXPAND(  Ex1_ave*Ex1_ave, 
                         + Ex2_ave*Ex2_ave, 
                         + Ex3_ave*Ex3_ave);
    #endif

    UU[k][j][i][ENG] += 0.5*(b2_new - b2_old);
    #endif 
  }
}
