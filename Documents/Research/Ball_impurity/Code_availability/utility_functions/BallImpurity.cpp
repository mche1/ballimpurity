#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double BallImpurity(int n, NumericMatrix D, IntegerVector I, double d)
  /* Input number of samples n, distance matrix D, intersting index of samples vector I and parameter d
   * Output the Ball Impurity of the subset I
   */
{
  int m=I.length(),i,j,k;
  int count=0,tempk;
  double BI=0;
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      count=0;
      for(k=0;k<m;k++)
      {
        tempk=I[k]-1;
        if(D(i,tempk) <= (D(i,j)+d))
        {count++;}
      }
      BI+=1.0*(1.0*count/m)*(1-1.0*count/m);
    }
  }
  BI=BI/(2*n*n);
  return(BI);
}
