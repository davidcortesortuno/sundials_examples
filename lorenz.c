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

/* NV_Ith_S: This macro gives access to the individual components of the
 * data array of an N Vector.  The assignment r = NV Ith S(v,i) sets r to
 * be the value of the i-th component of v.  Here i ranges from
 * 0 to n − 1 for a vector of length n */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
    /* Get the data from the y array: y = (F, W) */
    realtype theta = NV_Ith_S(y, 0);
    realtype omega = NV_Ith_S(y, 1);

    /* RHS of the equation for the derivative of W */
    realtype omegap = -sin(theta);

    /* Update the ydot vector: ydot = (D F, D W)
     *                                  t    t
     */
    NV_Ith_S(ydot,0) = omega;
    NV_Ith_S(ydot,1) = omegap;
    return 0;
}

int main(int argc, char** argv)
{
    int N = 200;
    realtype T0 = 0;
    realtype Tfinal = 100;
    realtype theta0 = atof(argv[1]);
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
    y = N_VNew_Serial(2);

    /* The assignment NV_Ith_S(v,i) = r sets the value of the i-th component of
     * v to be r.*/
    NV_Ith_S(y,0) = theta0; NV_Ith_S(y,1) = 0;

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

    /* In loop, call CVode, print results, and test for error. */
    for (k = 1; k < N; ++k) {
    realtype tout = k*Tfinal/N;
        if (CVode(cvode_mem, tout, y, &t, CV_NORMAL) < 0) {
            fprintf(stderr, "Error in CVode: %d\n", flag);
            return -1;
        }
        printf("%g %.16e %.16e\n", t, NV_Ith_S(y,0), NV_Ith_S(y,1));
    }
    N_VDestroy_Serial(y); /* Free y vector */
    CVodeFree(&cvode_mem); /* Free integrator memory */
    return(0);
}
