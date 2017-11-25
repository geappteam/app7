function [ h ] = identify( x, y, m )
% This function identifies a causal linear impulse response
% x is the input signal
% y is the output
% m is the size of the impulse response (nb Of coefficients)
% This functions returns m coefficients describing the impulse response

% Autocorrelation of original signal
Rxx = xcorr(x, x, m);

% Toepliz Matrix
A = toeplitz(Rxx(m+1:2*m+1));

% result table
Rxy = xcorr(x, y, m);
b = Rxy(1:m+1);

% Solving for h
h = linsolve(A, b);

end

