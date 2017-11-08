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

Note: also has the option to do (approximate) analytic priors.
All must be called with a string, and then following doubles (depending).
For now, four options:
  - Flat        : min, max (or, just max, assumes min=0)
  - SolidAngle  : takes no other arguments
  - Gaussian    : min, max
  - Log         : min, max
Note: these are all approximate. Still uses linear extrapolation between points.
*/

//******************************************************************************
NumericCdfInverse::NumericCdfInverse(std::string input_path_to_cdf)
/*
The constructor.
Takes in non-optional file, that holds numerical CDF, and forms the inverse.
See also overloaded version below.
*/
{

  //For the analytic "solid angle" (sin(theta)) prior
  if(input_path_to_cdf=="SolidAngle"){
    solidAnglePrior();
    return;
  }

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
NumericCdfInverse::NumericCdfInverse(std::string type, double a, double b)
/*
Overloaded initialisation function.
For (approximate) analytic priors.
Does Gaussian, and logarithmic priors.
*/
{

  //"solid angle" ALSO done above
  if(type=="SolidAngle"){
    solidAnglePrior();
  }
  else if(type=="Gaussian"){
    gaussianPrior(a,b);
  }
  else if(type=="log"||type=="Log"){
    logPrior(a,b);
  }
  else if(type=="flat"||type=="Flat"){
    flatPrior(a,b);
  }

}

//******************************************************************************
int NumericCdfInverse::flatPrior(double min, double max)
/*
Will lead to a flat prior between min/max.
Can just give 1 input, and will then go from 0 up to that number
*/
{
  if(max>min){
    xmin=min;
    xmax=max;
  }else{
    xmin=max;
    xmax=min;
  }
  inverse_cdf.push_back(xmin);
  inverse_cdf.push_back(xmax);
  N=2;
}

//******************************************************************************
int NumericCdfInverse::gaussianPrior(double x0, double s)
/*
Approximate Gaussian prior.
NOTE: not very accurate for very low/very high values for u..
but OK in most cases!
Only goes out to 4 sigma.
*/
{
  N=257; //must be an odd number!
  xmin=x0-4*s;
  xmax=x0+4*s;
  inverse_cdf.push_back(xmin);
  for(int i=1; i<N-1; i++){
    double u = double(i)/double(N-1);
    double g = x0 + s*sqrt(2)*inverseErf(2.*u-1.);
    inverse_cdf.push_back(g);
  }
  inverse_cdf.push_back(xmax);
  return 0;
}
//*-----------------------------------------------------------------------------
double NumericCdfInverse::inverseErf(double x)
/*
Approximate inverse Error function.
Note very accurate (~e-3).. but good enough for us here.
Needed for Gaussian prior.
*/
{
  int sgn=1;
  if(x<0)sgn=-1;
  double lnx=log(1.-x*x);
  double tt1=4.33+0.5*lnx;
  double tt2=6.803*lnx;
  return sgn*sqrt(sqrt(tt1*tt1-tt2)-tt1);
}

//******************************************************************************
int NumericCdfInverse::logPrior(double min, double max)
/*
Approximate log prior.
Note: still uses linear extrapolation between points..
good enough for most purposes though
*/
{
  N=256;
  xmin=fabs(min); //negative numbers not allowed
  xmax=fabs(max);
  for(int i=0; i<N; i++){
    double u = double(i)/double(N-1);
    double g = xmin*pow((xmax/xmin),u);
    inverse_cdf.push_back(g);
  }
  return 0;
}

//******************************************************************************
int NumericCdfInverse::solidAnglePrior()
/*
Approximate prior for the solid angle for the 'theta' angle
[z = cos(theta), theta: 0 -> Pi]
The prior is sin(theta)
*/
{
  N=128; //don't need many points.
  xmin=0;
  xmax=M_PI;
  for(int i=1; i<=N; i++){
    double u = double(i)/double(N);
    double g = acos(1. - 2.*u);
    inverse_cdf.push_back(g);
  }
  return 0;
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
