#ifndef spectralMaths_H
#define spectralMaths_H

/* This is a library for the spectral model k-distribution
It generates sets of k-values for a set of weights and g-values using
Gaussian quadrature. The number of quadrature points used is normally set
as 8.
There is also an interpolation routine provided

*/

#include <vector>

//function to generate Gaussian Quadrature and weights as per Chebyshave pol.
void quadgen(int, std::vector<double>, std::vector<double>);

//Quadrature generation as per Gauss Chebyshev scheme
void gausscheb2(std::vector<double>&,std::vector<double>&, int);

//generating k-values based on power law scheme with given power
std::vector<double> kPowerLaw(double,double, double, int);

// linear interpolation of k/g values from 64 to 8 points 
std::vector<double>
linearInterpMono
(
    int,
    const std::vector<double>&,
    const std::vector<double>&,
    int,
    const std::vector<double> &
);

#endif
