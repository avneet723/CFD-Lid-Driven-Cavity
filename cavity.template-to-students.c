/************************************************************************ */
/*      This code solves for the viscous flow in a lid-driven cavity      */
/**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/************* Following are fixed parameters for array sizes **************/
#define imax 9   	/* Number of points in the x-direction (use odd numbers only) */
#define jmax 9   	/* Number of points in the y-direction (use odd numbers only) */
#define neq 3       /* Number of equation to be solved ( = 3: mass, x-mtm, y-mtm) */

/**********************************************/
/****** All Global variables declared here. ***/
/***** These variables SHOULD not be changed **/
/********** by the program once set. **********/
/**********************************************/
/***** The variables declared "const" CAN *****/
/*** not be changed by the program once set ***/
/**********************************************/

/*--------- Numerical Constants --------*/
  const double zero   = 0.0;
  const double tenth  = 0.1;
  const double sixth  = 1.0/6.0;
  const double fifth  = 0.2;
  const double fourth = 0.25;
  const double third  = 1.0/3.0;
  const double half   = 0.5;
  const double one    = 1.0;
  const double two    = 2.0;
  const double three  = 3.0;
  const double four   = 4.0;
  const double six    = 6.0;
  
/*--------- User sets inputs here  --------*/

  const int nmax = 500000;        /* Maximum number of iterations */
  const int iterout = 5000;       /* Number of time steps between solution output */
  const int imms = 1;             /* Manufactured solution flag: = 1 for manuf. sol., = 0 otherwise */
  const int isgs = 1;             /* Symmetric Gauss-Seidel  flag: = 1 for SGS, = 0 for point Jacobi */
  const int irstr = 0;            /* Restart flag: = 1 for restart (file 'restart.in', = 0 for initial run */
  const int ipgorder = 0;         /* Order of pressure gradient: 0 = 2nd, 1 = 3rd (not needed) */
  const int lim = 0;              /* variable to be used as the limiter sensor (= 0 for pressure) */

  const double cfl  = 0.5;      /* CFL number used to determine time step */
  const double Cx = 0.01;     	/* Parameter for 4th order artificial viscosity in x */
  const double Cy = 0.01;      	/* Parameter for 4th order artificial viscosity in y */
  const double toler = 1.e-10; 	/* Tolerance for iterative residual convergence */
  const double rkappa = 0.1;   	/* Time derivative preconditioning constant */
  const double Re = 100.0;      	/* Reynolds number = rho*Uinf*L/rmu */
  const double pinf = 0.801333844662; /* Initial pressure (N/m^2) -> from MMS value at cavity center */
  const double uinf = 1.0;      /* Lid velocity (m/s) */
  const double rho = 1.0;       /* Density (kg/m^3) */
  const double xmin = 0.0;      /* Cavity dimensions...: minimum x location (m) */
  const double xmax = 0.05;   	/*                       maximum x location (m) */
  const double ymin = 0.0;      /*                       maximum y location (m) */
  const double ymax = 0.05;   	/*                       maximum y location (m) */
  const double Cx2 = 0.0;       /* Coefficient for 2nd order damping (not required) */
  const double Cy2 = 0.0;     	/* Coefficient for 2nd order damping (not required) */
  const double fsmall = 1.e-20; /* small parameter */

/*-- Derived input quantities (set by function 'set_derived_inputs' called from main)----*/
 
  double rhoinv =  -99.9; 	/* Inverse density, 1/rho (m^3/kg) */
  double rlength = -99.9;  	/* Characteristic length (m) [cavity width] */
  double rmu = -99.9;  		/* Viscosity (N*s/m^2) */
  double vel2ref = -99.9;  	/* Reference velocity squared (m^2/s^2) */
  double dx = -99.9; 		/*	 Delta x (m) */
  double dy = -99.9;  		/* Delta y (m) */
  double rpi = -99.9; 		/* Pi = 3.14159... (defined below) */

/*-- Constants for manufactured solutions ----*/
  const double phi0[neq] = {0.25, 0.3, 0.2};          /* MMS constant */
  const double phix[neq] = {0.5, 0.15, 1.0/6.0};      /* MMS amplitude constant */
  const double phiy[neq] = {0.4, 0.2, 0.25};          /* MMS amplitude constant */
  const double phixy[neq] = {1.0/3.0, 0.25, 0.1};     /* MMS amplitude constant */
  const double apx[neq] = {0.5, 1.0/3.0, 7.0/17.0}; 	/* MMS frequency constant */
  const double apy[neq] = {0.2, 0.25, 1.0/6.0};         /* MMS frequency constant */
  const double apxy[neq] = {2.0/7.0, 0.4, 1.0/3.0};     /* MMS frequency constant */
  const double fsinx[neq] = {0.0, 1.0, 0.0};            /* MMS constant to determine sine vs. cosine */
  const double fsiny[neq] = {1.0, 0.0, 0.0};            /* MMS constant to determine sine vs. cosine */
  const double fsinxy[neq] = {1.0, 1.0, 0.0};           /* MMS constant to determine sine vs. cosine */
                                                      /* Note: fsin = 1 means the sine function */
                                                      /* Note: fsin = 0 means the cosine function */
                                                      /* Note: arrays here refer to the 3 variables */ 
 
/********** Multidimensional array variables declared globally. *****************/
/******** Global declaration allows arrays to access larger amount of memory *******/
/***************** This enables running cases on larger meshes. ***************/
/***Use these variables cautiously as these are globally accessible from all functions.***/

  double u[imax][jmax][neq];         /* Solution vector [p, u, v]^T at each node */
  double uold[imax][jmax][neq];      /* Previous (old) solution vector */
  double s[imax][jmax][neq];    /* Source term */
  double dt[imax][jmax];        /* Local time step at each node */
  double artviscx[imax][jmax];  /* Artificial viscosity in x-direction */
  double artviscy[imax][jmax];  /* Artificial viscosity in y-direction */
  
/*------------- Function prototypes ----------------*/
/*----- All functions are globally accessible ------*/
  void set_derived_inputs();
  void output_file_headers();
  void initial(int* , double*, double resinit[neq]);
  void set_boundary_conditions();
  void bndry();
  void bndrymms();
  double umms(double, double, int);
  void write_output(int, double resinit[neq], double);
  void compute_source_terms();
  double srcmms_mass(double, double);  
  double srcmms_xmtm(double, double);
  double srcmms_ymtm(double, double);  
  void compute_time_step(double*);   
  void Compute_Artificial_Viscosity();
  void SGS_forward_sweep();
  void SGS_backward_sweep();
  void point_Jacobi();   
  void pressure_rescaling();    
  void check_iterative_convergence(int n, double res[neq], double resinit[neq], \
  int ninit, double rtime, double dtmin, double *conv);
  void Discretization_Error_Norms(double r1[neq], \
  	double r2[neq], double r3[neq]);  
  
/*--- Variables for file handling ---*/
/*--- All files are globally accessible ---*/
  
  FILE *fp1; /* For output of iterative residual history */
  FILE *fp2; /* For output of field data (solution) */
  FILE *fp3; /* For writing the restart file */
  FILE *fp4; /* For reading the restart file */  
  FILE *fp5; /* For output of final DE norms (only for MMS)*/  
//$$$$$$   FILE *fp6; /* For debug: Uncomment for debugging. */                      
/**************************************************************************/
/*      						Main Function					      	  */
/**************************************************************************/
int main()
{

/*----- Looping indices --------*/
  	
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  int k = 0;                       /* k index (# of equations) */
  int n = 0;	                   /* Iteration number index */

  double conv = -99.9 ; /* Minimum of iterative residual norms from three equations */

                                                      
/*--------- Solution variables declaration --------*/
  
  int ninit = 0;        	/* Initial iteration number (used for restart file) */

//$$$$$$   double u[imax][jmax][neq];         /* Solution vector [p, u, v]^T at each node */
//$$$$$$   double uold[imax][jmax][neq];      /* Previous (old) solution vector */
//$$$$$$   double s[imax][jmax][neq];    /* Source term */
//$$$$$$   double dt[imax][jmax];        /* Local time step at each node */
//$$$$$$   double artviscx[imax][jmax];  /* Artificial viscosity in x-direction */
//$$$$$$   double artviscy[imax][jmax];  /* Artificial viscosity in y-direction */
  double res[neq];              /* Iterative residual for each equation */
  double resinit[neq];          /* Initial iterative residual for each equation (from iteration 1) */
  double rL1norm[neq];          /* L1 norm of discretization error for each equation */     
  double rL2norm[neq];          /* L2 norm of discretization error for each equation */
  double rLinfnorm[neq];        /* Linfinity norm of discretization error for each equation */
  double rtime = -99.9;         /* Variable to estimate simulation time */
  double dtmin = 1.0e99;        /* Minimum time step for a given iteration (initialized large) */

  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */
  
/* Solution variables initialization with dummy values */
for (i=0; i<imax; i++)
{
 for (j=0; j<jmax; j++)
 {
     dt[i][j] = -99.9;
     artviscx[i][j] = -99.9;      
     artviscy[i][j] = -99.9;         
   for (k=0; k<neq; k++)
   {
     u[i][j][k] = -99.9;
     uold[i][j][k] = -99.9;
     s[i][j][k] = -99.9;
     res[k] = -99.9;
     resinit[k] = -99.9;
     res[k] = -99.9;
     rL1norm[k] = -99.9;
     rL2norm[k] = -99.9;
     rLinfnorm[k] = -99.9;        
   }
 }    
}

/* Debug output: Uncomment and modify if debugging */
//$$$$$$ fp6 = fopen("./Debug.dat","w");
//$$$$$$ fprintf(fp6,"TITLE = \"Debug Data Data\"\n");
//$$$$$$ fprintf(fp6,"variables=\"x(m)\"\"y(m)\"\"visc-x\"\"visc-y\"\n");
//$$$$$$ fprintf(fp6, "zone T=\"n=%d\"\n",n);
//$$$$$$ fprintf(fp6, "I= %d J= %d\n",imax, jmax);
//$$$$$$ fprintf(fp6, "DATAPACKING=POINT\n");

/* Set derived input quantities */
  set_derived_inputs();

/* Set up headers for output files */
  output_file_headers();

/* Set Initial Profile for u vector */
  initial(&ninit, &rtime, resinit);   

/* Set Boundary Conditions for u */
  set_boundary_conditions();

/* Write out inital conditions to solution file */
  write_output(ninit, resinit, rtime);
 
/* Initialize Artificial Viscosity arrays to zero (note: artviscx(i,j) and artviscy(i,j) */
for(i=0; i<imax; i++)
{
  for(j=0; j<jmax; j++)
  {
  	artviscx[i][j] = zero;
  	artviscy[i][j] = zero;
  }
}

/* Evaluate Source Terms Once at Beginning */
/*(only interior points; will be zero for standard cavity) */
  compute_source_terms();

/*========== Main Loop ==========*/
  for (n = ninit; n<= nmax; n++)
  {
    /* Calculate time step */  
	compute_time_step(&dtmin);

    /* Save u values at time level n (u and uold are 2D arrays) */
    for(i=0; i<imax; i++)
    {
      for(j=0; j<jmax; j++)
      {
        for(k=0; k<neq; k++)
        {
          uold[i][j][k] = u[i][j][k];
        }
      }
    }
    
    if(isgs==1) /* ==Symmetric Gauss Seidel== */
    {
      /* Artificial Viscosity */
      Compute_Artificial_Viscosity();
      
      /* Symmetric Gauss-Siedel: Forward Sweep */
      SGS_forward_sweep();
  
      /* Set Boundary Conditions for u */
      set_boundary_conditions();
   
      /* Artificial Viscosity */
      Compute_Artificial_Viscosity();
         	 
      /* Symmetric Gauss-Siedel: Backward Sweep */
      SGS_backward_sweep();

      /* Set Boundary Conditions for u */
      set_boundary_conditions();
    }
    else
    {
      if(isgs==0) /* ==Point Jacobi== */
      {
        /* Artificial Viscosity */
        Compute_Artificial_Viscosity();
      
      	/* Point Jacobi: Forward Sweep */
      	point_Jacobi();
   
      	/* Set Boundary Conditions for u */
      	set_boundary_conditions();
      }
      else
      {
      printf("ERROR: isgs must equal 0 or 1!\n");
      exit (0);     
      }
    }

    /* Pressure Rescaling (based on center point) */
    pressure_rescaling();

    /* Update the time */
    rtime = rtime + dtmin;

    /* Check iterative convergence using L2 norms of iterative residuals */
    check_iterative_convergence(n, res, resinit, ninit, \
    rtime, dtmin, &conv);

    if(conv<toler) 
    {
    fprintf(fp1, "%d %e %e %e %e\n",n, rtime, res[0], res[1], res[2]);
    goto converged;
    }
    
    /* Output solution and restart file every 'iterout' steps */
    if( ((n%iterout)==0) ) 
    {
     write_output(n, resinit, rtime);
    }
    
  }  /* ========== End Main Loop ========== */

  printf("Solution failed to converge in %d iterations!!!", nmax);
    
  goto notconverged;
    
converged:  /* go here once solution is converged */

  printf("Solution converged in %d iterations!!!", n);

notconverged:

  /* Calculate and Write Out Discretization Error Norms (will do this for MMS only) */
  Discretization_Error_Norms(rL1norm, rL2norm, rLinfnorm);

  /* Output solution and restart file */
  write_output(n, resinit, rtime);

  /* Close open files */
  fclose(fp1);
  fclose(fp2);
//$$$$$$   fclose(fp6); /* Uncomment for debug output */
}

/**************************************************************************/
/*      					All Other	Functions					      */
/**************************************************************************/
void set_derived_inputs()
{
  rhoinv = one/rho;                            /* Inverse density, 1/rho (m^3/kg) */
  rlength = xmax - xmin;                       /* Characteristic length (m) [cavity width] */
  rmu = rho*uinf*rlength/Re;                   /* Viscosity (N*s/m^2) */
  vel2ref = uinf*uinf;                         /* Reference velocity squared (m^2/s^2) */
  dx = (xmax - xmin)/(double)(imax - 1);          /* Delta x (m) */
  dy = (ymax - ymin)/(double)(jmax - 1);          /* Delta y (m) */
  rpi = acos(-one);                            /* Pi = 3.14159... */
  printf("rho,V,L,mu,Re: %f %f %f %f %f\n",rho,uinf,rlength,rmu,Re);
}
/**************************************************************************/
void output_file_headers()
{
  /*
  Uses global variable(s): imms, fp1, fp2
  */
  
  /* Note: The vector of primitive variables is: */
  /*               u = [p, u, v]^T               */  
  /* Set up output files (history and solution)  */    

    fp1 = fopen("./history.dat","w");
    fprintf(fp1,"TITLE = \"Cavity Iterative Residual History\"\n");
    fprintf(fp1,"variables=\"Iteration\"\"Time(s)\"\"Res1\"\"Res2\"\"Res3\"\n");

    fp2 = fopen("./cavity.dat","w");
    fprintf(fp2,"TITLE = \"Cavity Field Data\"\n");
    if(imms==1)
    {
      fprintf(fp2,"variables=\"x(m)\"\"y(m)\"\"p(N/m^2)\"\"u(m/s)\"\"v(m/s)\"");\
      fprintf(fp2,"\"p-exact\"\"u-exact\"\"v-exact\"\"DE-p\"\"DE-u\"\"DE-v\"\n");      
    }
    else
    {
      if(imms==0)
      {
        fprintf(fp2,"variables=\"x(m)\"\"y(m)\"\"p(N/m^2)\"\"u(m/s)\"\"v(m/s)\"\n");
      }      
      else
      {
        printf("ERROR! imms must equal 0 or 1!!!\n");
        exit (0);
      }       
    }

  /* Header for Screen Output */
  printf("Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum\n"); 
}
/**************************************************************************/
void initial(int* ninit, double* rtime, double resinit[neq])
{
  /* 
  Uses global variable(s): zero, one, irstr, imax, jmax, neq, uinf, pinf 
  To modify: ninit, rtime, resinit, u, s
   */
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  int k = 0;                       /* k index (# of equations) */
  
  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */

  /* This subroutine sets inital conditions in the cavity */

  /* Note: The vector of primitive variables is:  */
  /*              u = [p, u, v]^T               */
  
  if(irstr==0)   /* Starting run from scratch*/
  {  
    *ninit = 1;          /* set initial iteration to one */
    *rtime = zero;       /* set initial time to zero */
    for(k = 0; k<neq; k++)
    {
      resinit[k] = one;
    }
    for(i = 0; i<imax; i++)
    {
      for(j = 0; j<jmax; j++)
      {
        u[i][j][0] = pinf;
        u[i][j][1] = zero;
        u[i][j][2] = zero;
        s[i][j][0] = zero;
        s[i][j][1] = zero;
        s[i][j][2] = zero;
      }
    u[i][jmax-1][1] = uinf; /* Initialize lid (top) to freestream velocity */
    }
  }  
  else
  {
    if(irstr==1)  /* Restarting from previous run (file 'restart.in') */
    {
      fp4 = fopen("./restart.in","r"); /* Note: 'restart.in' must exist! */
	  if (fp4==NULL)
      {
		printf("Error opening restart file. Stopping.\n");
        exit (0);
      }      
      fscanf(fp4, "%d %lf", ninit, rtime); /* Need to known current iteration # and time value */
      fscanf(fp4, "%lf %lf %lf", &resinit[0], &resinit[1], &resinit[2]); /* Needs initial iterative residuals for scaling */
      for(j=0; j<jmax; j++)
      {
        for(i=0; i<imax; i++)
        {
         fscanf(fp4, "%lf %lf %lf %lf %lf", &x, &y, &u[i][j][0], &u[i][j][1], &u[i][j][2]); 
        }
      }
    *ninit = *ninit + 1;
    printf("Restarting at iteration %d\n", *ninit);
    fclose(fp4);
    }
    else
    {
      printf("ERROR: irstr must equal 0 or 1!\n");
      exit (0);
    }
  }
}
/**************************************************************************/
void set_boundary_conditions()
{
    /* 
  Uses global variable(s): imms
  To modify: u (via other functions: bndry() and bndrymms())
   */
     
  /* This subroutine determines the appropriate BC routines to call */
  if(imms==0) 
  {
    bndry();
  }
  else
  {
    if(imms==1)
    {
      bndrymms();
    }
    else
    {
      printf("ERROR: imms must equal 0 or 1!\n");
      exit (0);
    }
  }
}
/**************************************************************************/
void bndry()
{
    /* 
  Uses global variable(s): zero, one, two, half, imax, jmax, uinf  
  To modify: u 
   */
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */

  /* This applies the cavity boundary conditions */

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */


 
}
/**************************************************************************/
void bndrymms()
{
  /* 
  Uses global variable(s): two, imax, jmax, neq, xmax, xmin, ymax, ymin, rlength  
  To modify: u
  */
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  int k = 0;                       /* k index (# of equations) */
  
  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */

  /* This applies the cavity boundary conditions for the manufactured solution */

  /* Side Walls */
  for( j = 1; j<jmax-1; j++)
  {
    y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
    i = 0;
    x = xmin;
    for( k = 0; k<neq; k++)
    {
      u[i][j][k] = umms(x,y,k);
    }
    u[0][j][0] = two*u[1][j][0] - u[2][j][0];    /* 2nd Order BC */
//    u[0][j][0] = u[1][j][0];                  /* 1st Order BC */

    i=imax-1;
    x = xmax;
    for( k = 0; k<neq; k++)
    {
      u[i][j][k] = umms(x,y,k);
    }
    u[imax-1][j][0] = two*u[imax-2][j][0] - u[imax-3][j][0];   /* 2nd Order BC */
//	u[imax-1][j][0] = u[imax-2][j][0];                       /* 1st Order BC */
  }

/* Top/Bottom Walls */
 for(i=0; i<imax; i++)
 {
    x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
    j = 0;
    y = ymin;
    for( k = 0; k<neq; k++)
    {
      u[i][j][k] = umms(x,y,k);
    }
    u[i][0][0] = two*u[i][1][0] - u[i][2][0];   /* 2nd Order BC */
//$$$$$$     u[i][0][0] = u[i][1][0];            /* 1st Order BC */

    j = jmax-1;
    y = ymax;
    for( k = 0; k<neq; k++)
    {
      u[i][j][k] = umms(x,y,k);
    }
    u[i][jmax-1][0] = two*u[i][jmax-2][0] - u[i][jmax-3][0];   /* 2nd Order BC */
//$$$$$$     u[i][jmax-1][0] = u[i][jmax-2][0];              /* 1st Order BC */
 }
}
/**************************************************************************/
double umms(double x, double y, int k)  
{
  /* 
  Uses global variable(s): one, rpi, rlength
  Inputs: x, y, k
  To modify: <none>
  Returns: umms
  */

  double ummstmp; /* Define return value for umms as double precision */

  double termx = -99.9;      /* Temp variable */
  double termy = -99.9;      /* Temp variable */
  double termxy = -99.9;     /* Temp variable */
  double argx = -99.9;       /* Temp variable */
  double argy = -99.9;       /* Temp variable */
  double argxy = -99.9;      /* Temp variable */  
  
  /* This function returns the MMS exact solution */
  
  argx = apx[k]*rpi*x/rlength;
  argy = apy[k]*rpi*y/rlength;
  argxy = apxy[k]*rpi*x*y/rlength/rlength;
  termx = phix[k]*(fsinx[k]*sin(argx)+(one-fsinx[k])*cos(argx));
  termy = phiy[k]*(fsiny[k]*sin(argy)+(one-fsiny[k])*cos(argy));
  termxy = phixy[k]*(fsinxy[k]*sin(argxy)+(one-fsinxy[k])*cos(argxy));
  
  ummstmp = phi0[k] + termx + termy + termxy;
 
  return (ummstmp);  
}
/**************************************************************************/
void write_output(int n, double resinit[neq], double rtime)
{
   /* 
  Uses global variable(s): imax, jmax, new, xmax, xmin, ymax, ymin, rlength, imms
  Uses global variable(s): ninit, u, dt, resinit, rtime
  To modify: <none> 
  Writes output and restart files.
   */
   
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  int k = 0;                       /* k index (# of equations) */

  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */

  /* Field output */
  fprintf(fp2, "zone T=\"n=%d\"\n",n);
  fprintf(fp2, "I= %d J= %d\n",imax, jmax);
  fprintf(fp2, "DATAPACKING=POINT\n");

  if(imms==1) 
  {
    for(j=0; j<jmax; j++)
    {
      for(i=0; i<imax; i++)
      {
      x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
      y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
      fprintf(fp2,"%e %e %e %e %e %e %e %e %e %e %e\n", x, y, \
      u[i][j][0], u[i][j][1], u[i][j][2], umms(x,y,0), umms(x,y,1), umms(x,y,2), \
      (u[i][j][0]-umms(x,y,0)), (u[i][j][1]-umms(x,y,1)), (u[i][j][2]-umms(x,y,2)));
      }
    }    
  }
  else
  {
    if(imms==0)
    {
      for(j=0; j<jmax; j++)
      {
        for(i=0; i<imax; i++)
        {
        x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
        y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
        fprintf(fp2,"%e %e %e %e %e\n", x, y, \
        u[i][j][0], u[i][j][1], u[i][j][2]);
        }
      }
    }
    else
    {
      printf("ERROR: imms must equal 0 or 1!\n");
      exit (0);
    }
  }

  /* Restart file: overwrites every 'iterout' iteration */
  fp3 = fopen("./restart.out","w");    	
  fprintf(fp3,"%d %e\n", n, rtime);    
  fprintf(fp3,"%e %e %e\n", resinit[0], resinit[1], resinit[2]);
  for(j=0; j<jmax; j++)
  {
    for(i=0; i<imax; i++)
    {
    x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
    y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
    fprintf(fp3,"%e %e %e %e %e\n", x, y, \
    u[i][j][0], u[i][j][1], u[i][j][2]);
    }
  }
  fclose(fp3);
}
/**************************************************************************/
void compute_source_terms()
{
  /* 
  Uses global variable(s): imax, jmax, imms, rlength, xmax, xmin, ymax, ymin
  To modify: s (source terms)
  */

  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */

  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */

  /* Evaluate Source Terms Once at Beginning (only interior points; will be zero for standard cavity) */
  
  for(i=1; i<imax-1; i++)
  {
    for(j=1; j<jmax-1; j++)
    {
    x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
    y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
    s[i][j][0] = (double)(imms)*srcmms_mass(x,y);
    s[i][j][1] = (double)(imms)*srcmms_xmtm(x,y);
    s[i][j][2] = (double)(imms)*srcmms_ymtm(x,y);
    }
  }
}
/**************************************************************************/
double srcmms_mass(double x, double y)  
{
  /* 
  Uses global variable(s): rho, rpi, rlength
  Inputs: x, y
  To modify: <none>
  Returns: srcmms_mass
  */

  double srcmasstmp; /* Define return value for srcmms_mass as double precision */

  double dudx; 	/* Temp variable: u velocity gradient in x direction */
  double dvdy;  /* Temp variable: v velocity gradient in y direction */
  
/* This function returns the MMS mass source term */

  dudx = phix[1]*apx[1]*rpi/rlength*cos(apx[1]*rpi*x/rlength)  \
        + phixy[1]*apxy[1]*rpi*y/rlength/rlength  \
        * cos(apxy[1]*rpi*x*y/rlength/rlength);
  
  dvdy = -phiy[2]*apy[2]*rpi/rlength*sin(apy[2]*rpi*y/rlength)  \
       - phixy[2]*apxy[2]*rpi*x/rlength/rlength  \
       * sin(apxy[2]*rpi*x*y/rlength/rlength);

  srcmasstmp = rho*dudx + rho*dvdy;

  return (srcmasstmp);
} 
/**************************************************************************/
double srcmms_xmtm(double x, double y)  
{
  /* 
  Uses global variable(s): rho, rpi, rmu, rlength
  Inputs: x, y
  To modify: <none>
  Returns: srcmms_xmtm
  */

  double srcxmtmtmp; /* Define return value for srcmms_xmtm as double precision */

  double dudx; 	/* Temp variable: u velocity gradient in x direction */
  double dudy;  /* Temp variable: u velocity gradient in y direction */
  double termx;        /* Temp variable */
  double termy;        /* Temp variable */
  double termxy;       /* Temp variable */
  double uvel;         /* Temp variable: u velocity */
  double vvel;         /* Temp variable: v velocity */
  double dpdx;         /* Temp variable: pressure gradient in x direction */
  double d2udx2;       /* Temp variable: 2nd derivative of u velocity in x direction */
  double d2udy2;       /* Temp variable: 2nd derivative of u velocity in y direction */

/*This function returns the MMS x-momentum source term */

  termx = phix[1]*sin(apx[1]*rpi*x/rlength);
  termy = phiy[1]*cos(apy[1]*rpi*y/rlength);
  termxy = phixy[1]*sin(apxy[1]*rpi*x*y/rlength/rlength);
  uvel = phi0[1] + termx + termy + termxy;
  
  termx = phix[2]*cos(apx[2]*rpi*x/rlength);
  termy = phiy[2]*cos(apy[2]*rpi*y/rlength);
  termxy = phixy[2]*cos(apxy[2]*rpi*x*y/rlength/rlength);
  vvel = phi0[2] + termx + termy + termxy;
  
  dudx = phix[1]*apx[1]*rpi/rlength*cos(apx[1]*rpi*x/rlength) \ 
        + phixy[1]*apxy[1]*rpi*y/rlength/rlength  \
        * cos(apxy[1]*rpi*x*y/rlength/rlength);
  
  dudy = -phiy[1]*apy[1]*rpi/rlength*sin(apy[1]*rpi*y/rlength)  \
        + phixy[1]*apxy[1]*rpi*x/rlength/rlength  \
        * cos(apxy[1]*rpi*x*y/rlength/rlength);
  
  dpdx = -phix[0]*apx[0]*rpi/rlength*sin(apx[0]*rpi*x/rlength) \
        + phixy[0]*apxy[0]*rpi*y/rlength/rlength  \
        * cos(apxy[0]*rpi*x*y/rlength/rlength);

  d2udx2 = -phix[1]*pow((apx[1]*rpi/rlength),2)  \
          * sin(apx[1]*rpi*x/rlength)  \
          - phixy[1]*pow((apxy[1]*rpi*y/rlength/rlength),2)  \
          * sin(apxy[1]*rpi*x*y/rlength/rlength); 
 
  d2udy2 = -phiy[1]*pow((apy[1]*rpi/rlength),2)  \
          * cos(apy[1]*rpi*y/rlength)  \
          - phixy[1]*pow((apxy[1]*rpi*x/rlength/rlength),2)  \
          * sin(apxy[1]*rpi*x*y/rlength/rlength);
  
  srcxmtmtmp = rho*uvel*dudx + rho*vvel*dudy + dpdx  \
               - rmu*( d2udx2 + d2udy2 );

  return(srcxmtmtmp);
} 
/**************************************************************************/
double srcmms_ymtm(double x, double y)  
{
  /* 
  Uses global variable(s): rho, rpi, rmu, rlength
  Inputs: x, y
  To modify: <none>
  Returns: srcmms_ymtm
  */

  double srcymtmtmp; /* Define return value for srcmms_ymtm as double precision */

  double dvdx;         /* Temp variable: v velocity gradient in x direction */
  double dvdy;         /* Temp variable: v velocity gradient in y direction */
  double termx;        /* Temp variable */
  double termy;        /* Temp variable */
  double termxy;       /* Temp variable */
  double uvel;         /* Temp variable: u velocity */
  double vvel;         /* Temp variable: v velocity */
  double dpdy;         /* Temp variable: pressure gradient in y direction */
  double d2vdx2;       /* Temp variable: 2nd derivative of v velocity in x direction */
  double d2vdy2;       /* Temp variable: 2nd derivative of v velocity in y direction */

/* This function returns the MMS y-momentum source term */

  termx = phix[1]*sin(apx[1]*rpi*x/rlength);
  termy = phiy[1]*cos(apy[1]*rpi*y/rlength);
  termxy = phixy[1]*sin(apxy[1]*rpi*x*y/rlength/rlength);
  uvel = phi0[1] + termx + termy + termxy;
  
  termx = phix[2]*cos(apx[2]*rpi*x/rlength);
  termy = phiy[2]*cos(apy[2]*rpi*y/rlength);
  termxy = phixy[2]*cos(apxy[2]*rpi*x*y/rlength/rlength);
  vvel = phi0[2] + termx + termy + termxy;
  
  dvdx = -phix[2]*apx[2]*rpi/rlength*sin(apx[2]*rpi*x/rlength)  \
         - phixy[2]*apxy[2]*rpi*y/rlength/rlength  \
         * sin(apxy[2]*rpi*x*y/rlength/rlength);
  
  dvdy = -phiy[2]*apy[2]*rpi/rlength*sin(apy[2]*rpi*y/rlength)  \
         - phixy[2]*apxy[2]*rpi*x/rlength/rlength  \
         * sin(apxy[2]*rpi*x*y/rlength/rlength);
  
  dpdy = phiy[0]*apy[0]*rpi/rlength*cos(apy[0]*rpi*y/rlength)  \
         + phixy[0]*apxy[0]*rpi*x/rlength/rlength  \
         * cos(apxy[0]*rpi*x*y/rlength/rlength);
  
  d2vdx2 = -phix[2]*pow((apx[2]*rpi/rlength),2)  \
           * cos(apx[2]*rpi*x/rlength)  \
           - phixy[2]*pow((apxy[2]*rpi*y/rlength/rlength),2)  \
           * cos(apxy[2]*rpi*x*y/rlength/rlength);
  
  d2vdy2 = -phiy[2]*pow((apy[2]*rpi/rlength),2)  \
           * cos(apy[2]*rpi*y/rlength)  
           - phixy[2]*pow((apxy[2]*rpi*x/rlength/rlength),2)  \
           * cos(apxy[2]*rpi*x*y/rlength/rlength);
  
  srcymtmtmp = rho*uvel*dvdx + rho*vvel*dvdy + dpdy  \
                - rmu*( d2vdx2 + d2vdy2 );
  
  return (srcymtmtmp);  
}
/**************************************************************************/
void compute_time_step(double* dtmin)
{
  /* 
  Uses global variable(s): one, two, four, half, fourth
  Uses global variable(s): vel2ref, rmu, rho, dx, dy, cfl, rkappa, imax, jmax
  Uses: u
  To Modify: dt, dtmin
  */
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  
  double dtvisc = -99.9;      /* Viscous time step stability criteria (constant over domain) */
  double uvel2 = -99.9;       /* Local velocity squared */
  double beta2 = -99.9;       /* Beta squared paramete for time derivative preconditioning */
  double lambda_x = -99.9;    /* Max absolute value eigenvalue in (x,t) */
  double lambda_y = -99.9;    /* Max absolute value eigenvalue in (y,t) */
  double lambda_max = -99.9;  /* Max absolute value eigenvalue (used in convective time step computation) */
  double dtconv = -99.9;      /* Local convective time step restriction */

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */

  

}  
/**************************************************************************/
void Compute_Artificial_Viscosity()
{
  /* 
  Uses global variable(s): zero, one, two, four, six, half, fourth
  Uses global variable(s): imax, jmax, lim, rho, dx, dy, Cx, Cy, Cx2, Cy2, fsmall, vel2ref, rkappa
  Uses: u
  To Modify: artviscx, artviscy
  */
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */

  double uvel2 = -99.9;       /* Local velocity squared */
  double beta2 = -99.9;       /* Beta squared paramete for time derivative preconditioning */
  double lambda_x = -99.9;    /* Max absolute value e-value in (x,t) */
  double lambda_y = -99.9;    /* Max absolute value e-value in (y,t) */
  double d4pdx4 = -99.9;      /* 4th derivative of pressure w.r.t. x */
  double d4pdy4 = -99.9;      /* 4th derivative of pressure w.r.t. y */

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */

  
}
/**************************************************************************/
void SGS_forward_sweep()
{
  /* 
  Uses global variable(s): two, three, six, half
  Uses global variable(s): imax, imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, \
                        xmax, xmin, ymax, ymin, rmu, vel2ref
  Uses: artviscx, artviscy, dt, s
  To Modify: u
  */
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
 
  double dpdx = -99.9;        /* First derivative of pressure w.r.t. x */
  double dudx = -99.9;        /* First derivative of x velocity w.r.t. x */
  double dvdx = -99.9;        /* First derivative of y velocity w.r.t. x */
  double dpdy = -99.9;        /* First derivative of pressure w.r.t. y */
  double dudy = -99.9;        /* First derivative of x velocity w.r.t. y */
  double dvdy = -99.9;        /* First derivative of y velocity w.r.t. y */
  double d2udx2 = -99.9;      /* Second derivative of x velocity w.r.t. x */
  double d2vdx2 = -99.9;      /* Second derivative of y velocity w.r.t. x */
  double d2udy2 = -99.9;      /* Second derivative of x velocity w.r.t. y */
  double d2vdy2 = -99.9;      /* Second derivative of y velocity w.r.t. y */
  double beta2 = -99.9;       /* Beta squared parameter for time derivative preconditioning */
  double uvel2 = -99.9;       /* Velocity squared */

  /* Symmetric Gauss-Siedel: Forward Sweep */

  
/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




}
/**************************************************************************/
void SGS_backward_sweep()
{
  /* 
  Uses global variable(s): two, three, six, half
  Uses global variable(s): imax, imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, \
                        xmax, xmin, ymax, ymin, rmu, vel2ref
  Uses: artviscx, artviscy, dt, s
  To Modify: u
  */
  
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
 
  double dpdx = -99.9;        /* First derivative of pressure w.r.t. x */
  double dudx = -99.9;        /* First derivative of x velocity w.r.t. x */
  double dvdx = -99.9;        /* First derivative of y velocity w.r.t. x */
  double dpdy = -99.9;        /* First derivative of pressure w.r.t. y */
  double dudy = -99.9;        /* First derivative of x velocity w.r.t. y */
  double dvdy = -99.9;        /* First derivative of y velocity w.r.t. y */
  double d2udx2 = -99.9;      /* Second derivative of x velocity w.r.t. x */
  double d2vdx2 = -99.9;      /* Second derivative of y velocity w.r.t. x */
  double d2udy2 = -99.9;      /* Second derivative of x velocity w.r.t. y */
  double d2vdy2 = -99.9;      /* Second derivative of y velocity w.r.t. y */
  double beta2 = -99.9;       /* Beta squared parameter for time derivative preconditioning */
  double uvel2 = -99.9;       /* Velocity squared */

  /* Symmetric Gauss-Siedel: Backward Sweep  */
  
/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */



}
/**************************************************************************/
void point_Jacobi()
{
  /* 
  Uses global variable(s): two, three, six, half
  Uses global variable(s): imax, imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, \
                        xmax, xmin, ymax, ymin, rmu, vel2ref
  Uses: uold, artviscx, artviscy, dt, s
  To Modify: u
  */
  
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
 
  double dpdx = -99.9;        /* First derivative of pressure w.r.t. x */
  double dudx = -99.9;        /* First derivative of x velocity w.r.t. x */
  double dvdx = -99.9;        /* First derivative of y velocity w.r.t. x */
  double dpdy = -99.9;        /* First derivative of pressure w.r.t. y */
  double dudy = -99.9;        /* First derivative of x velocity w.r.t. y */
  double dvdy = -99.9;        /* First derivative of y velocity w.r.t. y */
  double d2udx2 = -99.9;      /* Second derivative of x velocity w.r.t. x */
  double d2vdx2 = -99.9;      /* Second derivative of y velocity w.r.t. x */
  double d2udy2 = -99.9;      /* Second derivative of x velocity w.r.t. y */
  double d2vdy2 = -99.9;      /* Second derivative of y velocity w.r.t. y */
  double beta2 = -99.9;       /* Beta squared parameter for time derivative preconditioning */
  double uvel2 = -99.9;       /* Velocity squared */

  /* Point Jacobi method */
  
/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */



}
/**************************************************************************/
void pressure_rescaling()
{
  /* 
  Uses global variable(s): imax, jmax, imms, xmax, xmin, ymax, ymin, rlength, pinf
  To Modify: u
  */
  
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */

  int iref = 0;                     /* i index location of pressure rescaling point */
  int jref = 0;                     /* j index location of pressure rescaling point */

  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */  
  double deltap = -99.9;  /* delta_pressure for rescaling all values */

  iref = (imax-1)/2;     /* Set reference pressure to center of cavity */
  jref = (jmax-1)/2;
  if(imms==1)
  {
    x = (xmax - xmin)*(double)(iref)/(double)(imax - 1);
    y = (ymax - ymin)*(double)(jref)/(double)(jmax - 1);
    deltap = u[iref][jref][0] - umms(x,y,0); /* Constant in MMS */
  }
  else
  {
    deltap = u[iref][jref][0] - pinf; /* Reference pressure */
  }

  for(i=0; i<imax; i++)
  {
    for(j=0; j<jmax; j++)
    {
    u[i][j][0] = u[i][j][0] - deltap;
    }
  }      
}  
/**************************************************************************/
void check_iterative_convergence(int n, double res[neq], double resinit[neq], \
int ninit, double rtime, double dtmin, double *conv)
{
  /* 
  Uses global variable(s): zero
  Uses global variable(s): imax, jmax, neq, fsmall
  Uses: n, u, uold, dt, res, resinit, ninit, rtime, dtmin
  To modify: conv
  */
   
  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  int k = 0;                       /* k index (# of equations) */

  /* Compute iterative residuals to monitor iterative convergence */
  
/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




  /* Write iterative residuals every 10 iterations */
  if( ((n%10)==0)||(n==ninit) )
  {
    fprintf(fp1, "%d %e %e %e %e\n",n, rtime, res[0], res[1], res[2] );
    printf("%d   %e   %e   %e   %e   %e\n",n, rtime, dtmin, res[0], res[1], res[2] );    
    /* Maybe a need to format this better */    
  }
      
  /* Write header for iterative residuals every 200 iterations */
  if( ((n%200)==0)||(n==ninit) )
  {
  printf("Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum\n"); 
  }  
}
/**************************************************************************/
void Discretization_Error_Norms(double rL1norm[neq], \
double rL2norm[neq], double rLinfnorm[neq]) 
{
  /* 
  Uses global variable(s): zero
  Uses global variable(s): imax, jmax, neq, imms, xmax, xmin, ymax, ymin, rlength
  Uses: u
  To modify: rL1norm, rL2norm, rLinfnorm 
  */

  int i = 0;                       /* i index (x direction) */
  int j = 0;                       /* j index (y direction) */
  int k = 0;                       /* k index (# of equations) */

  double x = -99.9;       /* Temporary variable for x location */
  double y = -99.9;       /* Temporary variable for y location */
  double DE = -99.9;  	/* Discretization error (absolute value) */

  if(imms==1)
  {
  
/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */



  }
  
}
