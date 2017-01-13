/*
* Lorenz equation:
*
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

/* System variables */
#define SIGMA 10
#define RHO 28
#define BETA 8/3

/* Type : UserData (contains grid constants) */
typedef struct {
  realtype sigma, rho, beta;
} *UserData;

/* Taken from one of the CVODE examples */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;
  
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
      return(1); }
  
    /* Check if flag < 0 */
    else if (opt == 1) {
      errflag = (int *) flagvalue;
      if (*errflag < 0) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                funcname, *errflag);
        return(1); }}
  
    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
      fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
      return(1); }
  
    return(0);
}

/* NV_Ith_S: This macro gives access to the individual components of the
 * data array of an N Vector.  The assignment r = NV Ith S(v,i) sets r to
 * be the value of the i-th component of v.  Here i ranges from
 * 0 to n − 1 for a vector of length n */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    UserData data;

    /* Extract needed constants from data */
    data = (UserData) user_data;
    realtype S = data->sigma;
    realtype R = data->rho;
    realtype B = data->beta;

    /* Get the data from the y array: y = (F, W) */
    realtype X = NV_Ith_S(y, 0);
    realtype Y = NV_Ith_S(y, 1);
    realtype Z = NV_Ith_S(y, 2);

    /* Update the ydot vector: ydot = ( D X, D Y, D Z )
     *                                   t    t    t
     */
    NV_Ith_S(ydot, 0) = S * (Y - X);
    NV_Ith_S(ydot, 1) = X * (R - Z) - Y;
    NV_Ith_S(ydot, 2) = X * Y - B * Z;

    return 0;
}

int main(int argc, char** argv)
{
    UserData data;

    /* Number of points and final-initial times */
    int N = 5000;
    realtype T0 = 0;
    realtype Tfinal = 50;

    /* Integrator tolerances */
    realtype reltol = 1e-6;
    realtype abstol = 1e-8;
    
    realtype t;
    int flag, k;
    N_Vector y = NULL;
    void* cvode_mem = NULL;

    /* Create serial vector of length NEQ for I.C.
     * Manual: This function creates and allocates memory for a serial N
     * Vector. Its only argument is the vector length.
     * */
    y = N_VNew_Serial(3);

    /* INITIAL VALUES
     * The assignment NV_Ith_S(v,i) = r sets the value of the i-th component of
     * v to be r.*/
    NV_Ith_S(y, 0) = 1;
    NV_Ith_S(y, 1) = 2;
    NV_Ith_S(y, 2) = 1;

    /* LORENZ EQS COEFFICIENTS */
    data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
    if(check_flag((void *)data, "malloc", 2)) return(1); 
    data->sigma = SIGMA;
    data->rho   = RHO;
    data->beta  = BETA;

    /* Set up solver */
    cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    if (cvode_mem == 0) {
        fprintf(stderr, "Error in CVodeMalloc: could not allocate\n");
        return -1;
    }

    /* Call CVodeMalloc to initialize the integrator memory */
    /* flag = CVodeMalloc(cvode_mem, f, T0, y, CV_SS, reltol, &abstol); */
    flag = CVodeInit(cvode_mem, f, T0, y);
    if (flag < 0) {
        fprintf(stderr, "Error in CVodeMalloc: %d\n", flag);
        return -1;
    }

    flag = CVodeSStolerances(cvode_mem, reltol, abstol);

    /* Set the pointer to user-defined data */
    flag = CVodeSetUserData(cvode_mem, data);
    if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

    /* In loop, call CVode, print results, and test for error. */
    for (k = 1; k < N; ++k) {
        /* Time: */
        realtype tout = k * Tfinal / N;
        /* Run until time */
        if (CVode(cvode_mem, tout, y, &t, CV_NORMAL) < 0) {
            fprintf(stderr, "Error in CVode: %d\n", flag);
            return -1;
        }
        /* Print results (better to define a file output) */
        printf("%g %.16e %.16e %.16e\n", t, NV_Ith_S(y,0),
                                            NV_Ith_S(y,1),
                                            NV_Ith_S(y,2));
    }

    N_VDestroy_Serial(y); /* Free y vector */
    CVodeFree(&cvode_mem); /* Free integrator memory */
    return(0);
}
