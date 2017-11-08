//class PriorCdfInverse::
#include "NumericCdfInverseClass.h"

/*
171005.
Contains member functions for the class NumericCdfInverse

Class that finds the inverse of numerical CDFs.
Used for inverse transform sampling, when the inverse either cannot be found
analytically or is too complicated.

This class reads in a given numerical CDF function.
It then forms the inverse of CDF, by solving numerically.
It takes advantage of the fact that the CDF is a smooth and monotonically
increasing function; therefore this 'solve/invert' routine definitely won't
work for any general function.

Then, there is a single public function, inverseCdf, that outputs the value of
the inverse CDF function for any given u=[0,1].
It uses linear extrapolation.


*/

//******************************************************************************
NumericCdfInverse::NumericCdfInverse(std::string input_path_to_cdf)
/*
The constructor.
Takes in non-optional file, that holds numerical CDF, and forms the inverse.

NOTE:
Could make a second optional input, an N-dimentional double vector, that holds
the coeficients of an (N-1)-dimentional polynomial.
i.e., the 'best-fit' polynomial for any arbitrary
NOTE: this would ONLY work, if the polynomial was guarenteed to be monotonically
increasing... prob best not to bother.

*/
{
  //Store the path the the CDF file:
  path_to_cdf = input_path_to_cdf;
  //sanity check:
  ok=false;
  //Read in the CDF file:
  readNumericCdf();
  //Use CDF to generate g(u) [the inverse CDF]
  invertCdf();
  if(N<=0){
    std::cout<<"\n\n!XX!\n FAIL 46 in NumericCdfInverse: "<<input_path_to_cdf
             <<" - No CDF supplied??\n!XX!\n\n";
  }else{
    ok=true;
  }
}



//******************************************************************************
NumericCdfInverse::NumericCdfInverse(double min, double max)
/**/
{

  xmin=min;
  xmax=max;
  inverse_cdf.push_back(xmin);
  inverse_cdf.push_back(xmax);
  N=2;

}



//******************************************************************************
int NumericCdfInverse::readNumericCdf()
/*
Open and read in the numeric CDF [Cumulative Distribution Function] for the
probability distribution
Note: CDF function must be in a plain text file, no headers or comments.
Format should be each line has 'x' 'space' 'f(x)', e.g.
  x_1 y_1
  x_2 y_2
  ...
  x_N y_n
The First x (x_min) is typically 0, but not always.
The last x is x_max, and depends on the function.
The first y should always be 0, and the last y should always be 1.
Note: Each proceeding y MUST be larger than previous y!

*/
{

  //Open the input file that containts the CDF:
  std::ifstream cdf_file;
  cdf_file.open (path_to_cdf.c_str());

  //Read through each line of file, write into array
  std::string temp_line;
  bool first_line = true;
  while ( getline(cdf_file,temp_line) ) {
    std::stringstream ss_line (temp_line);
    double x,y;
    ss_line >> x >> y;
    if(first_line){
      // Find xmin (the first x given)
      xmin=x; //just the first time
      first_line=false;
    }
    xmax=x; // Will continuously update, until the end!
    //nb: for now, only works for LINEAR spacing! Can update for non-uniform!
    //Store value in array:
    cdf.push_back(y);
  }

  //get number of points, and step-size:
  N = (int) cdf.size();
  dx = (xmax - xmin) / (N - 1);

  cdf_file.close();

  return 0;
}


//******************************************************************************
int NumericCdfInverse::invertCdf()
/*
Finds the inverse of the given CDF.
Stores it in the array inverse_cdf.
Solves the equation
  cdf(x) = u
for u = [0,1]
Uses a linear extrapolation.
*/
{

  double ixm = 0;  // x_minus (integer index for cdf array)
  //nb: x_minus = ixm * dx + xmin

  //Loop through each value of u, find x such that cdf(x)=u
  for(int i=0; i<N; i++){
    double u = double(i)/(N-1); // u: o->1

    // Solve the equation cdf(x)=u for x
    // Since cdf is smooth and monotonically increasing, we don't need to start
    // from scratch each time!
    int ixp = ixm; // x_plus; integer index for x
    while( cdf[ixp] <= u ){
      ixp++;
      if(ixp>=N){
        //this should never happen..
        ixp = N-1;
        break;
      }
    }
    if(ixp!=0) ixm = ixp - 1;

    //use linear extrapolation between x_minus and x_plus
    double a = u - cdf[ixm];
    double b = cdf[ixp] - u;
    double x = xmin; // if a=0 and b=0 [happens for u=0], use x_min (?)
    if(a+b>0)
      x = xmin + dx*(b*ixm + a*ixp) / (a+b);

    //for safety (I think this never happens?)
    if(x>xmax)x=xmax;
    if(x<xmin)x=xmin;

    //Store in the inverse cdf array:
    inverse_cdf.push_back(x);

  }

  return 0;
}



//******************************************************************************
double NumericCdfInverse::inverseCdf(double u)
/*
This is the public function.
It returns:
  x = g(u) = inverse_cdf(u)
Note: inverse_cdf is an array.
Uses linear extrapolation between two nearest integers.
*/
{
  // short-cut for min/max values:
  if(u<=0) return xmin;
  if(u>=1) return xmax;

  // Map u = [0,1]  to  du = [0,N-1]
  double diu = (N-1)*u;

  // Find nearest two integer values for array index
  int ium = int(diu); //i_u_minus
  int iup = ium + 1;  //i_u_plus
  if(iup >= N ) iup = N-1; //don't overshoot array! (delta=0 in this case!)
  double delta = diu - ium;

  // Output the x=g(u) value, using linear extrapolation:
  return inverse_cdf[ium] * (1-delta) + inverse_cdf[iup] * delta;
}
