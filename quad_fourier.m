function J = quad_fourier(g,CS,a,b,w,n,d)
%QUAD_FOURIER Compute the integral of the function
%       g(x)*x^(a-1)*exp(-b*x)*cos(w*x)    (1)
% or
%       g(x)*x^(a-1)*exp(-b*x)*sin(w*x)    (2)
% over [0, Inf), by using the coupled Gaussian rule developed in [1].
%   J = QUAD_FOURIER(G,CS,A,B,W,N) returns in J the approximate value over
%   [0,Inf) of the integral of (1) or (2) with parameters A,B,W, by using
%   N quadrature points. CS = 'cos' is relative to (1), while CS = 'sin' 
%   refers to function (2).
%   J = QUAD_FOURIER(G,CS,A,B,W,N,D) allows the user to work with variable
%   precision arithmetic with D signifincant decimal digits. The output
%   is returned in standard double precision.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
%   [1] E. Denich and P. Novati, "Numerical quadrature for integrals
%       involving oscillating functions", 2024.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the inputs
if nargin < 6
    error('Not enough input arguments')
end

if nargin == 6
    d = [];
end

f = @(t) g(t/w); c = b/w;

[xF,wF] = fourierpts(n,a,c,CS,d); 
[xL,wL] = lagpts(n,a-1); 

IF = sum(f(xF).*wF');
IL = 1/c^a*sum(f(xL./c).*wL');

I = IF-IL;
J = I/w^a;

end