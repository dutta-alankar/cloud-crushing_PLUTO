/* ----------- Auto generated from generateCatalystAdaptor.py ----------------- 
Created on Sat Mar 30 23:55:52 2024

@author: alankar                                                                */

#include <catalyst.h>
#include <stdio.h>
#include <math.h>
#include "pluto.h"

//-----------------------------------------------------------------------------
/**
 * Initialize Catalyst.
 */
//-----------------------------------------------------------------------------
void do_catalyst_initialization(int scripts, char** pipeline_script)
{
  conduit_node* catalyst_init_params = conduit_node_create();
  // pass scripts pass on command line.
  char buf[1024];
  int i;
  for (i=0; i<scripts; i++){
    snprintf(buf, 256, "catalyst/scripts/script%d", i);
    conduit_node_set_path_char8_str(catalyst_init_params, buf, pipeline_script[i]);
  }
  conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/implementation", "paraview");
  conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/search_paths/paraview", PARAVIEW_IMPL_DIR);
  enum catalyst_status err = catalyst_initialize(catalyst_init_params);
  conduit_node_destroy(catalyst_init_params);
  if (err != catalyst_status_ok)
  {
    printf("Failed to initialize Catalyst: %d\n", err);
  }
}

//-----------------------------------------------------------------------------
/**
 * Execute per cycle
 */
//-----------------------------------------------------------------------------
void do_catalyst_execute(int cycle, double time, Grid* grid, Data* d)
{
  /*   Conversion to Conduit (catalyst) compatible data arrangement   */
  
  unsigned int i, j, k;
  static unsigned int initialized = 0;
  unsigned int counter = 0;
  
  static double *xp, *yp, *zp;
  static uint64_t *CellConn;
  const uint64_t numCells = grid->np_tot[IDIR]*grid->np_tot[JDIR]*grid->np_tot[KDIR];
  const uint64_t numPoints[3] = {grid->np_tot[IDIR]+1, grid->np_tot[JDIR]+1, grid->np_tot[KDIR]+1};
  
  //double ***temperature  = GetUserVar("temperature");
  //ComputeUserVar (d, grid);  

  if (initialized==0){
    xp = (double *)malloc((grid->np_tot[IDIR]+1) * (grid->np_tot[JDIR]+1) * (grid->np_tot[KDIR]+1) * sizeof(double));
    yp = (double *)malloc((grid->np_tot[IDIR]+1) * (grid->np_tot[JDIR]+1) * (grid->np_tot[KDIR]+1) * sizeof(double));
    zp = (double *)malloc((grid->np_tot[IDIR]+1) * (grid->np_tot[JDIR]+1) * (grid->np_tot[KDIR]+1) * sizeof(double));
    CellConn   = (uint64_t*)malloc(8 * numCells * sizeof(int64_t));
    
    for(k=0; k<grid->np_tot[KDIR]; k++)
    {
      for(j=0; j<grid->np_tot[JDIR]; j++)
      {
        for(i=0; i<grid->np_tot[IDIR]; i++)
        { 
          //cell_count = k * grid->np_tot[IDIR] * grid->np_tot[JDIR] + j * grid->np_tot[IDIR] + i;
          CellConn[counter++] =     k * numPoints[JDIR]*numPoints[IDIR] +     j * numPoints[IDIR] + i; 

          CellConn[counter++] = (k+1) * numPoints[JDIR]*numPoints[IDIR] +     j * numPoints[IDIR] + i;

          CellConn[counter++] = (k+1) * numPoints[JDIR]*numPoints[IDIR] + (j+1) * numPoints[IDIR] + i;

          CellConn[counter++] =     k * numPoints[JDIR]*numPoints[IDIR] + (j+1) * numPoints[IDIR] + i;;

          CellConn[counter++] =     k * numPoints[JDIR]*numPoints[IDIR] +     j * numPoints[IDIR] + i+1;

          CellConn[counter++] = (k+1) * numPoints[JDIR]*numPoints[IDIR] +     j * numPoints[IDIR] + i+1;

          CellConn[counter++] = (k+1) * numPoints[JDIR]*numPoints[IDIR] + (j+1) * numPoints[IDIR] + i+1;

          CellConn[counter++] =     k * numPoints[JDIR]*numPoints[IDIR] + (j+1) * numPoints[IDIR] + i+1;
        }
      }
    }
  }
  
  counter = 0;
  double x1, x2, x3;
  for (k = 0; k<=grid->np_tot[KDIR]; k++){
    for (j = 0; j<=grid->np_tot[JDIR]; j++){
      for (i = 0; i<=grid->np_tot[IDIR]; i++){
        x1 = (i!=grid->np_tot[IDIR])?grid->xl[IDIR][i]:grid->xr[IDIR][i-1];
        x2 = (j!=grid->np_tot[JDIR])?grid->xl[JDIR][j]:grid->xr[JDIR][j-1];
        x3 = (k!=grid->np_tot[KDIR])?grid->xl[KDIR][k]:grid->xr[KDIR][k-1];

	#if GEOMETRY==CARTESIAN || GEOMETRY == CYLINDRICAL
	DIM_EXPAND(
        xp[counter]   = x1;,
        yp[counter]   = x2;,
        zp[counter]   = x3;
        )
	#elif GEOMETRY==SPHERICAL
	  #if DIMENSIONS <= 2
	DIM_EXPAND(
	xp[counter] = x1*sin(x2);,
	yp[counter] = x1*cos(x2);,
	zp[counter] = 0.;
	)
	  #else
	DIM_EXPAND(
        xp[counter]   = x1*sin(x2)*cos(x3);,
        yp[counter]   = x1*sin(x2)*sin(x3);,
        zp[counter]   = x1*cos(x2);
        )
          #endif
        #elif GEOMETRY==POLAR
	DIM_EXPAND(
        xp[counter]   = x1*cos(x2);,
        yp[counter]   = x1*sin(x2);,
        zp[counter]   = x3;
        )
        #else
        printLog ("! CatalystAdaptor: Unknown geometry\n");
        QUIT_PLUTO(1);
	#endif
	counter++;
      }
    }
  }
  counter = 0;

  /*   setup conduit for catalyst   */
  conduit_node* catalyst_exec_params = conduit_node_create();
  conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep", cycle);
  conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/cycle", cycle);
  conduit_node_set_path_float64(catalyst_exec_params, "catalyst/state/time", time);
  
  conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh");
  conduit_node* mesh = conduit_node_create();
  
  // add coordsets
  conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "explicit");
  
  conduit_node_set_path_external_float64_ptr(mesh, "coordsets/coords/values/X", xp, numPoints[IDIR]*numPoints[JDIR]*numPoints[KDIR]);
  conduit_node_set_path_external_float64_ptr(mesh, "coordsets/coords/values/Y", yp, numPoints[IDIR]*numPoints[JDIR]*numPoints[KDIR]);
  conduit_node_set_path_external_float64_ptr(mesh, "coordsets/coords/values/Z", zp, numPoints[IDIR]*numPoints[JDIR]*numPoints[KDIR]);
  
  // add topologies
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "unstructured");
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords");
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/elements/shape", "hex");
  conduit_node_set_path_external_int64_ptr(mesh, "topologies/mesh/elements/connectivity", CellConn, numCells * 8);

  // add rho (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/rho/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/rho/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/rho/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/rho/values", /*rho */(double *)(**d->Vc[RHO]), numCells );    
    
  // add prs (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/prs/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/prs/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/prs/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/prs/values", /*prs */(double *)(**d->Vc[PRS]), numCells );    
    
  // add vx1 (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/vx1/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/vx1/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/vx1/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/vx1/values", /*vx1 */(double *)(**d->Vc[VX1]), numCells );    
    
  // add vx2 (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/vx2/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/vx2/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/vx2/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/vx2/values", /*vx2 */(double *)(**d->Vc[VX2]), numCells );    
    
  // add vx3 (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/vx3/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/vx3/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/vx3/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/vx3/values", /*vx3 */(double *)(**d->Vc[VX3]), numCells );    
    
  // add tr1 (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/tr1/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/tr1/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/tr1/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/tr1/values", /*tr1 */(double *)(**d->Vc[TRC]), numCells );    
    
  ComputeUserVar (d, grid);
    
  // add Temp (cell-field)
  double ***Temp  = GetUserVar("Temp");
  conduit_node_set_path_char8_str(mesh, "fields/Temp/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/Temp/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/Temp/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/Temp/values", /*Temp */(double *)(**Temp), numCells ); 
        
  // add ndens (cell-field)
  double ***ndens  = GetUserVar("ndens");
  conduit_node_set_path_char8_str(mesh, "fields/ndens/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/ndens/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/ndens/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/ndens/values", /*ndens */(double *)(**ndens), numCells ); 
        
  // add mach (cell-field)
  double ***mach  = GetUserVar("mach");
  conduit_node_set_path_char8_str(mesh, "fields/mach/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/mach/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/mach/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/mach/values", /*mach */(double *)(**mach), numCells ); 
        
  // add PbykB (cell-field)
  double ***PbykB  = GetUserVar("PbykB");
  conduit_node_set_path_char8_str(mesh, "fields/PbykB/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/PbykB/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/PbykB/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/PbykB/values", /*PbykB */(double *)(**PbykB), numCells ); 
        
  // add cellvol (cell-field)
  double ***cellvol  = GetUserVar("cellvol");
  conduit_node_set_path_char8_str(mesh, "fields/cellvol/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/cellvol/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/cellvol/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr(mesh, "fields/cellvol/values", /*cellvol */(double *)(**cellvol), numCells ); 
        
  // add the mesh info (conduit mesh) to catalyst_exec_params
  conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh);
  
  #ifdef CATALYST_DEBUG
  // print for debugging purposes, if needed
  conduit_node_print(catalyst_exec_params);

  // print information with details about memory allocation
  conduit_node* info = conduit_node_create();
  conduit_node_info(catalyst_exec_params, info);
  conduit_node_print(info);
  conduit_node_destroy(info);
  #endif
  
  initialized = 1;
  enum catalyst_status err = catalyst_execute(catalyst_exec_params);
  if (err != catalyst_status_ok)
  {
    printf("Failed to execute Catalyst: %d\n", err);
  }
  conduit_node_destroy(catalyst_exec_params);
  conduit_node_destroy(mesh);
}

//-----------------------------------------------------------------------------
/**
 * Finalize Catalyst.
 */
//-----------------------------------------------------------------------------
void do_catalyst_finalization()
{
  conduit_node* catalyst_fini_params = conduit_node_create();
  enum catalyst_status err = catalyst_finalize(catalyst_fini_params);
  if (err != catalyst_status_ok)
  {
    printf("Failed to execute Catalyst: %d\n", err);
  }
  conduit_node_destroy(catalyst_fini_params);
}
