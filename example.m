close all
clear al

% Example of usage of quad_fourier.m
g = @(x) 1./(1+exp(-x)); CS = 'cos';
a = 0.5; b = 0.4; w = 2; n = 40; d = 30;

Je = quad_fourier_lag(g,CS,a,b,w,200); % reference solution
J = quad_fourier(g,CS,a,b,w,n,d);

disp('Relative error:  ')
disp(abs((J-Je)/Je))
