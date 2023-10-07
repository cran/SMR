/* https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#C   */

#include <Rcpp.h>
#include <vector>

double Pi;

using namespace std;
using namespace Rcpp;

void lege_coef(NumericMatrix *lcoef, int N)
{
  int n, i;
  (*lcoef)(0, 0) = (*lcoef)(1, 1) = 1;
  for (n = 2; n <= N; n++)
  {
    (*lcoef)(n, 0) = -(n - 1) * (*lcoef)(n - 2, 0) / n;
    for (i = 1; i <= n; i++)
      (*lcoef)(n, i) = ((2 * n - 1) * (*lcoef)(n - 1, i - 1) - (n - 1) * (*lcoef)(n - 2, i)) / n;
  }
}

double lege_eval(NumericMatrix *lcoef, int n, double x)
{
  int i;
  double s = (*lcoef)(n, n);
  for (i = n; i; i--)
    s = s * x + (*lcoef)(n, i - 1);
  return s;
}

double lege_diff(NumericMatrix *lcoef, int n, double x)
{
  return n * (x * lege_eval(lcoef, n, x) - lege_eval(lcoef, n - 1, x)) / (x * x - 1);
}

void lege_roots(std::vector<double> &lroots, NumericMatrix *lcoef, std::vector<double> &weight, int N)
{
  int i;
  double x, x1;
  for (i = 1; i <= N; i++)
  {
    x = cos(Pi * (i - .25) / (N + .5));
    do
    {
      x1 = x;
      x -= lege_eval(lcoef, N, x) / lege_diff(lcoef, N, x);
    } while (fdim(x, x1) > 2e-16);
    /*  fdim( ) was introduced in C99, if it isn't available
     *  on your system, try fabs( ) */
    lroots[i - 1] = x;

    x1 = lege_diff(lcoef, N, x);
    weight[i - 1] = 2 / ((1 - x * x) * x1 * x1);
  }
}

double lege_inte(double *weight, int N, double *lroots, double (*f)(double), double a, double b)
{
  double c1 = (b - a) / 2, c2 = (b + a) / 2, sum = 0;
  int i;
  for (i = 0; i < N; i++)
    sum += weight[i] * f(c1 * lroots[i] + c2);
  return c1 * sum;
}
//[[Rcpp::export]]
NumericMatrix gaussLegendre(int n)
{

  Pi = atan2(1, 1) * 4;
  int i, j;
  std::vector<double> lroots(n, 0);
  NumericMatrix out(2, n);
  std::vector<double> weight(n, 0);
  // double lcoef[n + 1][n + 1] = {{0}};
  NumericMatrix lcoef(n + 1);
  //	double ** lcoef = new double * [1 + n];
  //	for (i = 0; i < 1 + n; i++) {
  //		lcoef[i] = new double[1 + n];
  //	}
  for (i = 0; i < 1 + n; i++)
  {
    for (j = 0; j < 1 + n; j++)
    {
      lcoef(i, j) = 0;
    }
  }

  lege_coef(&lcoef, n);
  lege_roots(lroots, &lcoef, weight, n);

  for (i = 0; i < n; i++)
  {
    out(0, i) = lroots[i];
    out(1, i) = weight[i];
  }

  return out;
}
