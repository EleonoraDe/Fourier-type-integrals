function [x,w] = fourierpts(n,a,c,CS,d)
%FOURIERPTS Nodes and weights relative to the weight function
%
%             w_C(t) = t^(a-1)*exp(-c*t)*(cos(t)+1),
% or
%             w_S(t) = t^(a-1)*exp(-c*t)*(sin(t)+1).
%
%   [X,W] = FOURIERPTS(N,A,C,CS) returns N nodes in (0,Inf) in the column
%   vector X and a row vector W of weights, relative to the weight function
%   w_C(t) or w_S(t), with parameters A and C. A and C must be real-valued
%   scalar > 0. CS = 'cos' is relative to the weight function w_C(t),
%   while CS = 'sin' allows the user to work with w_S(t).
%
%   [X,W] = FOURIERPTS(N,A,C,CS,D) allows the user to work with variable
%   precision arithmetic with D significant decimal digits. The outputs
%   are returned in standard double precision.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTES AND REFERENCES:
%
% Methods:
% The code is adapted from [1].
% The Chebyshev algorithm is adapted from [2].
% The Golub-Welsch algorithm is adapted from [3].
%
% References:
%   [1] E. Denich and P. Novati, "Numerical quadrature for integrals
%       involving oscillating functions", 2024.
%   [2] W. Gautschi, "On generating orthogonal polynomials", SIAM J. SCl.
%       Stat. Comput., 3(3):289-317, 1969.
%   [3] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the inputs

if n < 1 
    error('First input should be a positive integer.')
end

if  a <= 0
    error('Second input should be positive.')
end

if c <= 0
    error('Third input should be positive.')
end

if ~ strcmp(CS,'cos') && ~ strcmp(CS,'sin')
    error('Fourth input should be cos or sin.')
end

if nargin == 5
    syms mu 
    a = sym(a); c = sym(c);
else
    d = [];
end

% Computation of the moments 
phi = atan(1/c);
if CS == 'cos'
    F = cos(a*phi);
else
    F = sin(a*phi);
end
s = sqrt(1+c^2);
mu(1) = gamma(a)/s^a*F+gamma(a)/c^a;
for k = 1:2*n-1
    if CS == 'cos'
        G = cos((k+a)*phi);
        H = cos((k-1+a)*phi);
    else
        G = sin((k+a)*phi);
        H = sin((k-1+a)*phi);
    end
    mu(k+1) = (k-1+a)/c*(G*(cos(phi))^(k+a)+1)/...
        (H*(cos(phi))^(k-1+a)+1)*mu(k);
end

% Computation of the recurrence coefficients via Chebyshev algorithm
ab=chebyshev(n,mu,d);

% Computation of nodes and weights using Golub-Welsch algorithm
mu1 = double(mu(1));
alpha=ab(1:n,1);
beta=ab(1:n,2);
rbeta=sqrt(beta); 
J=diag(alpha)+diag(rbeta(2:end),-1)+diag(rbeta(2:end),1);
[V,D] = eig(J);                       
x = diag(D);                           
[x, i] = sort(x);                     
w = mu1*V(1,i).^2;                     

% Routine for the Chebyshev algorithm

function ab=chebyshev(n,mu,d)

if ~isempty(d)
    syms sig
    digits(d);
    mu=vpa(mu,d);
end

ab(1,1)=mu(2)/mu(1); ab(1,2)=mu(1);
sig(1,1:2*n)=0; sig(2,:)=mu(1:2*n);
for N=3:n+1
  for m=N-1:2*n-N+2
    sig(N,m)=sig(N-1,m+1)-ab(N-2,1)*sig(N-1,m)-ab(N-2,2)*sig(N-2,m);
  end
  ab(N-1,1)=sig(N,N)/sig(N,N-1)-sig(N-1,N-1)/sig(N-1,N-2);
  ab(N-1,2)=sig(N,N-1)/sig(N-1,N-2);
end
ab=double(ab);
end

end
