#include "NumericCdfInverseClass.h"

int main()
{

  NumericCdfInverse cdf_object("testCdf.txt");  //numeric
  //NumericCdfInverse cdf_object(3,5);  // 3 -> 5
  //NumericCdfInverse cdf_object(200);    // 0 -> 200

  for( double u=0; u<=1; u+=(1./16) ){
    // Outputting in steps of '1/16' to demonstrate correctness
    // even when the extrapolation is used
    std::cout<<u<<" "<<cdf_object.inverseCdf(u)<<"\n";
  }

}
