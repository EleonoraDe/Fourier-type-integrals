This little package contains the main routine fourierpts.m for the computation of nodes and weights relative to the weigth function
w(t) = t^(alpha-1) * exp(-c*t) * (cos(t)+1) or w(t) = t^(alpha-1) * exp(-c*t) * (sin(t)+1), in [0, +Inf). 
The function quad_fourier.m evaluates the integral of the function g(x)*x^(alpha-1)*exp(-beta*x)*cos(w*x) or g(x)*x^(alpha-1)*exp(-beta*x)*sin(w*x), by
using the coupled Gaussian rule developed in [1].

The function quad_fourier_lag.m evaluates the integral of the same functions, by using the generalized Gauss-Laguerre rule.
Both codes require the function lagpts.m of the Chebfun open-source software [2].

An example is shown in the script example.m, in which the reference solution is computed by using the code quad_fourier_lag.m.
All the codes are implemented in Matlab.


[1] E. Denich and P. Novati, "Numerical quadrature for integal involving oscillating functions", 2024.

[2] T. A. Driscoll, N. Hale, and L. N. Trefethen, editors, Chebfun Guide, Pafnuty Publications, Oxford, 2014.
