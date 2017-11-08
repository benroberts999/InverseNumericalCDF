#ifndef _INVCDF_H
#define _INVCDF_H
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>

//Class that finds the inverse of numerical CDFs.
//Used for inverse transform sampling, when the inverse either cannot be
// found analytically or is too complicated.

class NumericCdfInverse {

  public:

    //constructor:
    NumericCdfInverse(std::string input_path_to_cdf);
    NumericCdfInverse(std::string type, double a, double b=0);

    //Function to return the inverse CDF value:
    double inverseCdf(double x);

    //for sanity check:
    bool ok;

  private:

    //Function that reads in the input CDF file:
    int readNumericCdf();

    // Function that inverts the given CDF numerically
    int invertCdf();

    //Functions for the "analytic" priors
    int solidAnglePrior();
    int flatPrior(double in_min, double in_max);
    int gaussianPrior(double x0, double s);
    double inverseErf(double x);
    int logPrior(double min, double max);

    //Store location of CDF input file:
    std::string path_to_cdf;

    // CDF domain:
    double xmin;  //u=0 shhould map to this
    double xmax;  //u=1 should map to this
    int N;        // Number of steps in numeric CDF
    double dx;    // CDF domain step size

    //Vectors that store the cdf and its inverse:
    std::vector<double> cdf;
    std::vector<double> inverse_cdf;

};

#endif
