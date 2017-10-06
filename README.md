# Inverse Numerical CDF

Class that finds the inverse of numerical CDFs (cumulative distribution function).

The resultant inverse CDF can then be used for inverse transform sampling
  see: https://github.com/benroberts999/inverseTransformSample

The reason you might want to form the inverse from the CDF numerically, is that
it is not always possible to find the inverse analytically.
And even in many cases where it is possible, it's often not practical!

By definition, a cumulative distribution function has range:

  CDF(x) = [0,1]

for

  x = [x_min, m_max]

Then define g(u) to be the inverse of CDF(x)  [solve CDF(x) = u]

  g(u) = [x_min, x_max]

for 

  u = [0,1].

**Inverse Transform Sampling Reminder:** 
If u is drawn randomly from [0,1] (with a uniform distribution), and CDF(x) is the CDF for probability distribution P(x), then g(u) will have probability distribution P(x), where g is the inverse of CDF.
Therefore, if we want to draw random numbers according to P(x), we just draw uniform random numbers u=[0,1], and then use g(u).


### What the program does:

This class reads in a given numerical CDF function (from a text file).
It then forms the inverse of the CDF, by solving numerically.
It takes advantage of the fact that the CDF is a smooth and monotonically
increasing function; therefore this 'solve/invert' routine definitely won't
work for a general function.

Then, there is a single public function, inverseCdf, that outputs the value of
the inverse CDF function for any given u=[0,1].
It uses linear extrapolation, which is almost always good enough.

The more points that are given in the input CDF file, the more accurate the output inverse CDF will be.

The CDF function must be in a plain text file, no headers or comments.
Format should be each line has 'x' 'space' 'f(x)', e.g.  
  x_1 y_1  
  x_2 y_2  
  ...  
  x_N y_n  

The First x (x_min) is typically 0, but not always.
The last x is x_max, and depends on the function.
Each proceeding y MUST be larger than previous y!

### An example:

The "exampleProgram.cpp" just shows a very basic example.

An example numerical CDF is given in the text file "testCdf.txt"
This CDF is for a Gaussian with mean 10, and standard deviation 3.
The CDF is given for x: 0 -> 20, which is enough.

Note: in this case, the inverse can actually be found analytically, allowing
one to easily check the correctness:

 g(u) = 10 + 3 Sqrt(2) InverseErf( 2 * u - Erf[ 5 * Sqrt{2} ] )
