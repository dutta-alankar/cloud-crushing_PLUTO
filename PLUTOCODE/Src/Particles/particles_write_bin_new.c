/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Writer for particles binary data in .dbl, .flt and 
        ASCII in .tab (only for serial version) format .
 
 \authors A. Mignone (mignone@to.infn.it)\n
          B. Vaidya (bvaidya@unito.it)\n
 
 \date    Aug 17, 2020
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void Particles_WriteHeader (PHeader, Output *, unsigned long int);

/* ********************************************************************************** */
void Particles_WriteBinary(particleNode *PHeadRef, double dt_particles, Output *output)
/*! 
 *  Write particle data in single or double precision binary format.
 *  The binary file structure is:
 *
 *  <Header (ASCII) section>
 *   .
 *   .
 *   .
 *  {field1, field2, ...}_p
 *   .
 *   .
 *   . 
 *
 * All fields are mapped from the particle structure into an array
 * and thus converted to either single or double precision, depending
 * on the function that is being called.
 * Fields can be scalars or array: for each particle nelem indicates
 * the number of elements written for each particle.
 * The field order does not have to match the structure member order but it
 * must be consistent with the sequence declared in Particles_SetOutput()
 * function.
 * Individual fields may be excluded by calling SetOutputVar() 
 * from ChangeOutputVar ().
 *
 *  \param [in]  PheadRef      Pointer to the Head Node of Particle List.
 *  \param [in]  dt_particles  Particle time step
 *  \param [in]  output        Pointer to output structure
 *
 ************************************************************************* */
{
  size_t   size;
  int      dir, nv, k, nfields, nelem;
  int     *dump_var = output->dump_var;
  int      nvar      = output->nvar;
  long int i, offset, offset_bytes, nparticles_glob, proc_npart[g_nprocs+1];
  void    *arr;
  float   *farr;
  double  *darr, u[3], gamma;
  unsigned long int npmax, npartWrite, leftToWrite, cnt;
  int nfile = 0, l, TotFiles, np_g, nf, np;
  char tempchar[128], filename[128];
 
  PHeader Header;
  particleNode * CurNode;

  
/* --------------------------------------------------------
   1. Loop over particles and map struct. members to array
   -------------------------------------------------------- */
  if (output->type == PARTICLES_DBL_OUTPUT)
      npmax = MAXNPART_PER_FILE_DBL;

  if (output->type == PARTICLES_FLT_OUTPUT)
      npmax = MAXNPART_PER_FILE_FLT;

  offset      = 0;
  leftToWrite = 0;
#ifdef PARALLEL
  MPI_Allreduce(&p_nparticles, &nparticles_glob, 1,
                MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(&p_nparticles, 1, MPI_LONG, &(proc_npart[1]), 1,
                MPI_LONG, MPI_COMM_WORLD);

  /* Compute individual processor offset (in particle units) */
  for(l = 0; l < prank; l++) {
      if (offset + proc_npart[l+1] > npmax) {
          leftToWrite = offset + proc_npart[l+1] - npmax;
          nfile++;
          offset = 0;
          while (leftToWrite > npmax){/* remaining particles over many files */
                 nfile++;
                 leftToWrite -= npmax ;
               }
          offset = leftToWrite; /*set offset to remainder */
        }else offset += proc_npart[l+1];
      }
#else
  nparticles_glob = p_nparticles;
#endif

  nfields = 0; /* Count how many fields are written to disk */
  nelem   = 0; /* The number of elements to be written (some fields may
                  be arrays). */
  for (nv = 0; nv < nvar; nv++) {
      if (dump_var[nv]) {
         nfields++;
         nelem += output->field_dim[nv];
        }}

 /* ------------------------------------------
  *  1.       Write the header 
  *  -----------------------------------------*/
  Header.npart_g  = nparticles_glob;
  Header.dt_p     = dt_particles;
  Header.nfields  = nfields;

  Particles_WriteHeader (Header, output, npmax);
  

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  np_g  = nparticles_glob;
  TotFiles  = np_g/npmax + 1;
  if ((np_g % npmax == 0) && (np_g > 0)) TotFiles -= 1;

  CurNode     = PHeadRef;
  leftToWrite = p_nparticles;

 /*---------------------------------------------
  * 2.    Loop over TotFiles
  * --------------------------------------------*/ 
  for (nf = 0; nf < TotFiles; nf++){

  sprintf (filename,"%s/particles.%04d", output->dir, output->nfile);
  sprintf(tempchar,"_%02d.%s",nf,output->ext);
  strcat(filename,tempchar);

  printLog ("> Writing particle file #%d (%s) to disk...", 
          output->nfile, output->ext);

/* Open files */
  if (output->type == PARTICLES_DBL_OUTPUT) size = sizeof(double);
  if (output->type == PARTICLES_FLT_OUTPUT) size = sizeof(float);
#ifdef PARALLEL
  MPI_File fhw;
  MPI_Status status;
  MPI_Offset start;   /* Store end of file position */
  
  MPI_File_open(MPI_COMM_WORLD, filename, 
               MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fhw);
  MPI_File_get_position(fhw, &start);

  offset_bytes = offset*nelem*size; /* Convert offset to bytes */
  MPI_File_set_view (fhw, start+offset_bytes, MPI_BYTE, MPI_CHAR, "native", MPI_INFO_NULL);
#else  
  FILE *fp = fopen(filename,"ab");
#endif 

  if (nf == nfile){

 /* ------------------------------------------
  *  3.      Begin loop on local particles 
  *  -----------------------------------------*/
  
 
      if (offset + leftToWrite > npmax){
          npartWrite   = npmax - offset;
          leftToWrite -= npartWrite;
          /* Set nfile & offset for next file */
          nfile++;
          offset = 0;
         }  
      else npartWrite = leftToWrite; 

     darr = ARRAY_1D(nelem*npartWrite, double);
     if (output->type == PARTICLES_FLT_OUTPUT) 
        farr = ARRAY_1D(nelem*npartWrite, float);

     i = 0; /* Initialize array index */

     for (np = 0; np < npartWrite; np++){

     /* -- Compute three- or four-velocity -- */
        for (dir = 0; dir < 3; dir++) u[dir] = CurNode->p.speed[dir];

  /* ------------------------------------------------------
         Map structure members to array.
         Important: field order should match the order 
         given in Particles_SetOutput().
         Here nv scan *all* fields (nv <= nfields)
     ------------------------------------------------------ */

       nv = 0;
       #if (PARTICLES == PARTICLES_CR) || (PARTICLES == PARTICLES_DUST)
       if (dump_var[nv++]) darr[i++] = CurNode->p.id;
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
       if (dump_var[nv++]) darr[i++] = u[IDIR];
       if (dump_var[nv++]) darr[i++] = u[JDIR];
       if (dump_var[nv++]) darr[i++] = u[KDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.mass;
       if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
       if (dump_var[nv++]) darr[i++] = CurNode->p.color;
       #endif

       #if PARTICLES == PARTICLES_DUST
       if (dump_var[nv++]) darr[i++] = CurNode->p.id;
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
       if (dump_var[nv++]) darr[i++] = u[IDIR];
       if (dump_var[nv++]) darr[i++] = u[JDIR];
       if (dump_var[nv++]) darr[i++] = u[KDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.mass;
       if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
       if (dump_var[nv++]) darr[i++] = CurNode->p.color;
       #endif

       #if PARTICLES == PARTICLES_LP
       if (dump_var[nv++]) darr[i++] = CurNode->p.id;
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
       if (dump_var[nv++]) darr[i++] = u[IDIR];
       if (dump_var[nv++]) darr[i++] = u[JDIR];
       if (dump_var[nv++]) darr[i++] = u[KDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
       if (dump_var[nv++]) {
           for (k = 0; k < PARTICLES_NCOLORS; k++) darr[i++] = CurNode->p.color[k];
           }
       #if PARTICLES_LP_SPECTRA == YES
       if (dump_var[nv++]) darr[i++] = CurNode->p.density;  // Density is the fluid
                                                         // density interpolated
                                                         // at the particle position
                                                         // and it is only needed to get N from chi
       if (dump_var[nv++]) darr[i++] = CurNode->p.pressure;  
       if (dump_var[nv++]) darr[i++] = CurNode->p.mag[IDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.mag[JDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.mag[KDIR];
       if (dump_var[nv++]) darr[i++] = CurNode->p.nmicro;
       if (dump_var[nv++]) darr[i++] = CurNode->p.cmp_ratio; // Could be eliminated from writing
                                                          // And also from structure (can be computed
                                                          // one at a time).
       if (dump_var[nv++]) darr[i++] = CurNode->p.shkflag;   
       if (dump_var[nv++]) darr[i++] = CurNode->p.shk_gradp;
       if (dump_var[nv++]) darr[i++] = CurNode->p.ca;        // Could be eliminated from writing since
       if (dump_var[nv++]) darr[i++] = CurNode->p.cr;        // they're computed from fluid quantities

       if (dump_var[nv++]) {
          for (k = 0; k < NFLX; k++) darr[i++] = CurNode->p.Vshk_upst[k];
          }

       if (dump_var[nv++]) {
          for (k = 0; k < NFLX; k++) darr[i++] = CurNode->p.Vshk_dnst[k];
          }

       if (dump_var[nv++]){
          for(k = 0; k <= PARTICLES_LP_NEBINS; k++) {
             darr[i++] = CurNode->p.eng[k];
         }}
    
       if (dump_var[nv++]){
         for(k = 0; k < PARTICLES_LP_NEBINS; k++){
           darr[i++] = CurNode->p.chi[k];
         }}

       #endif  /* PARTICLES_LP_SPECTRA == YES  */
       #endif  /* PARTICLES == PARTICLES_LP */
 
      CurNode = CurNode->next;
     }
  

/* --------------------------------------------------------
   4. Write data 
   -------------------------------------------------------- */
 
     if (output->type == PARTICLES_DBL_OUTPUT) size = sizeof(double);
     if (output->type == PARTICLES_FLT_OUTPUT) size = sizeof(float);

     if (output->type == PARTICLES_FLT_OUTPUT){
        for (i = 0; i < nelem*npartWrite; i++) farr[i] = (float)darr[i];
        arr = (void *) farr;
      }else if (output->type == PARTICLES_DBL_OUTPUT)
              arr = (void *) darr;

     char *cbuf = (char *) arr;
     cnt = nelem*npartWrite*size;
     #ifdef PARALLEL
        MPI_File_write(fhw, cbuf, cnt, MPI_CHAR, &status);
     #else 
        fwrite (cbuf, sizeof(char), cnt, fp);
     #endif 

/* --------------------------------------------------------
   5. Free_Array 
   -------------------------------------------------------- */

    FreeArray1D(darr);   
    if (output->type == PARTICLES_FLT_OUTPUT)
        FreeArray1D(farr);

  }//--nf==nfile
  
  #ifdef PARALLEL
     MPI_File_close(&fhw);
     MPI_Barrier(MPI_COMM_WORLD);
  #else
     fclose (fp);
  #endif

 }//--nfile
 
}


/***************************************************************************** */
void Particles_WriteHeader (PHeader Header, Output *output, unsigned long int npmax)
/* 
 * Write file header for all files
 *
 * *************************************************************************** */
{
  char fheader[1024], tempchar[128], filename[128];
  int  *dump_var = output->dump_var;
  int  nvar      = output->nvar;
  int  i, nf, npart=0, TotFiles;
  unsigned long int np_g;

  np_g = Header.npart_g;
  TotFiles    = np_g/npmax + 1;
  if ((np_g % npmax == 0) && np_g > 0) TotFiles -= 1;

  for (nf = 0; nf < TotFiles; nf++){
      sprintf (filename,"%s/particles.%04d", output->dir, output->nfile);
      sprintf(tempchar,"_%02d.%s",nf,output->ext);
      strcat(filename,tempchar);

      npart = MIN(npmax, np_g - nf*npmax);       
 
      sprintf(fheader,"# PLUTO %s binary particle data file    \n", PLUTO_VERSION);
      sprintf(fheader+strlen(fheader),"# dimensions     %d     \n", DIMENSIONS);
      sprintf(fheader+strlen(fheader),"# nflux          %d     \n", NFLX);
      sprintf(fheader+strlen(fheader),"# dt_particles   %12.6e \n", Header.dt_p);  
      if (IsLittleEndian())
          sprintf(fheader+strlen(fheader),"# endianity      little\n");
      else
          sprintf(fheader+strlen(fheader),"# endianity      big\n");

      sprintf(fheader+strlen(fheader),"# nparticles     %d     \n", npart);  
      sprintf(fheader+strlen(fheader),"# Totparticles   %ld    \n", np_g);  
      sprintf(fheader+strlen(fheader),"# nfiles         %d     \n", TotFiles);  
      sprintf(fheader+strlen(fheader),"# idCounter      %ld    \n", p_idCounter);
      sprintf(fheader+strlen(fheader),"# particletype   %d     \n", PARTICLES);
      if (output->type == PARTICLES_FLT_OUTPUT){  
          sprintf(fheader+strlen(fheader),"# precision      float\n");    
      }else if (output->type == PARTICLES_DBL_OUTPUT){  
          sprintf(fheader+strlen(fheader),"# precision      double\n");    
      }
  
      sprintf(fheader+strlen(fheader),"# time           %12.6e \n", g_time);  
      sprintf(fheader+strlen(fheader),"# stepNumber     %ld    \n", g_stepNumber);
      sprintf(fheader+strlen(fheader),"# nfields        %d     \n", Header.nfields);
      sprintf(fheader+strlen(fheader),"# field_names    ");
      for (i = 0; i < nvar; i++){
          if (dump_var[i]) {
              sprintf(fheader+strlen(fheader),output->var_name[i]);
              sprintf(fheader+strlen(fheader),"  ");
             }  
          }
      sprintf(fheader+strlen(fheader),"\n");

      sprintf(fheader+strlen(fheader),"# field_dim      ");
      for (i = 0; i < nvar; i++){
          if (dump_var[i]) {
             sprintf(fheader+strlen(fheader),"%d",output->field_dim[i]);
             sprintf(fheader+strlen(fheader),"  ");
          }}  
  
      sprintf(fheader+strlen(fheader),"\n");

      FileWriteHeader(fheader, filename, -1); /* NOTE: FileWriteHeader mode -1, 
                                                       creates new file */

    }
}


/* ********************************************************************* */
void Particles_WriteTab(particleNode* PHeadRef, char filename[128])
/*
 * Write particle coordinates, ids and speeds in Tab ASCII format 
 * only for *serial* version. 
 *
 *  \param [in]  PheadRef    Pointer to the Head Node of Particle List.
 *  \param [in]  filename    File name of particle data: particles.nnnn.tab,
 *                           where nnnn is the data file number.
 *********************************************************************** */
{
#ifdef PARALLEL
  printLog("! WARNING: Particle Data in Tabulated ASCII format is only written \
            with serial version");
#else
  int i;
  FILE *stream;
  particleNode* CurNode;
  
  stream = fopen(filename, "w");
  fprintf(stream, "# Nparticles: %ld\n", p_nparticles);
  fprintf(stream, "# Step:       %ld\n", g_stepNumber);
  fprintf(stream, "# Time:       %f\n",  g_time);
  
  
  CurNode = PHeadRef;

  while(CurNode != NULL) {
    fprintf(stream, "%d", CurNode->p.id);
    
    for (i = 0; i < 3; ++i) {
      fprintf(stream, "  %lf", CurNode->p.coord[i]);
    }
    for (i = 0; i < 3; ++i) {
      fprintf(stream, "  %lf", CurNode->p.speed[i]);
    }
    fprintf(stream, "\n");

    CurNode = CurNode->next;
  }
  fclose(stream);
#endif
}

