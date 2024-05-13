function J = quad_fourier_lag(g,CS,a,b,w,n)
%QUAD_FOURIER_LAG Compute the integral of the function
%       g(x)*x^(a-1)*exp(-b*x)*cos(w*x)    (1)
% or
%       g(x)*x^(a-1)*exp(-b*x)*sin(w*x)    (2)
% over [0, Inf), by using the generalized Gauss-Laguerre rule.
%   J = QUAD_FOURIER_LAG (G,CS,A,B,W,N) returns in J the approximate value 
%   over [0,Inf) of the integral of (1) or (2) with parameters A,B,W, by 
%   using N quadrature points. CS = 'cos' is relative to (1), while 
%   CS = 'sin' refers to function (2).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTES AND REFERENCES:
% The nodes and weights of the generalized Gauss-Laguerre rule are computed
% by using the function lagpts.m of the Chebfun package (see [2]).
%
% References:
%   [1] E. Denich and P. Novati, "Numerical quadrature for integrals
%       involving oscillating functions", 2024.
%   [2] T. A. Driscoll, N. Hale, and L. N. Trefethen, editors, Chebfun 
%       Guide, Pafnuty Publications, Oxford, 2014.

% Check the inputs
if nargin < 6
    error('Not enough input arguments')
end

if CS == 'cos'
    f = @(y) g(y/b).*cos(w/b*y);
else
    f = @(y) g(y/b).*sin(w/b*y);
end
[xL,wL] = lagpts(n,a-1); 

J = 1/b^a*sum(f(xL).*wL');


end