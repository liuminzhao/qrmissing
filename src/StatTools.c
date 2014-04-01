#include <R.h>
#include <Rmath.h>

/* 2012/09/13 add rbeta */

double F77_SUB(invcdfnorm)(double *p, double *mu, double *sigma, int *lower_tail, int *log_p)
{
  return qnorm(*p, *mu, *sigma, *lower_tail, *log_p);
}

double F77_SUB(dnrm)(double *x, double *mu, double *sigma, int *give_log)
{
  return dnorm(*x, *mu, *sigma, *give_log);
}

double F77_SUB(pnrm)(double *x, double *mu, double *sigma, int *lower_tail, int *give_log)
{
  return pnorm(*x, *mu, *sigma, *lower_tail, *give_log);
}


double F77_SUB(dgamma2)(double *x, double *shape, double *scale, int *give_log)
{
  return dgamma(*x, *shape, *scale, *give_log);
}


void F77_SUB(rndstart)(void) { GetRNGstate(); }

void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(myrnorm)(double *mu, double *sigma) { return rnorm(*mu, *sigma) ;}
double F77_SUB(myrunif)(double *a, double *b) { return runif(*a, *b) ;}

double F77_SUB(myrbeta)(double *a, double *b) { return rbeta(*a, *b) ; }

double F77_SUB(mydbeta)(double *x, double *a, double *b, int *give_log){
  return dbeta(*x, *a, *b, *give_log);
}

int F77_SUB(ihmssf)(int *i, int *j, int *n){
  int ans;
  int i1;
  int j1;

  i1=*i-1;
  j1=*j-1;

  if (i <= j)
    {
      ans=*n*i1-*i*i1/2+*j ;
    }
  else
    {
      ans=*n*j1-*j*j1/2+*i;
    }

  return ans;
}

double F77_SUB(mylddirichlet)(double *x, double *alpha, int *k){
  /* give log pdf of dirichlet distribution of x with parameter alpha and k as the length */
  int i;
  double ans;

  ans = 0;
  for (i = 0; i < *k; i ++){
    ans += (alpha[i] - 1) * log(x[i]);
  }

  return ans;
}
