function output = SumESF(X, N, list)
%-=========================================================================
% Implementation of the sum of elementary symmetric functions.
%
% Let  X = ( x_1 , . . . , x_M ) and let S(X) = sum{i1<...<iN}
% x_i1*x_i2*<...*x_iN such that N<=M. Then this function computes S(X),
% also known as the sum of elementary symmetric functions.
%
% The algorithm that carries out the sum is taken from: “Accurate validated
% and fast evaluation of elementary symmetric functions and its
% application,” Hao Jiang et al, 2016.
%
% Inputs:------------------------------------------------------------------
%   X: is a Mx1 dimentional array as explained above.
%
%   N: is the number of products in each summand. Must be an integer.
%
%   list: Contains the indices of those elements that need to be eliminated
%   from x. One dimensional array of integers.
%--------------------------------------------------------------------------
% Output:
%   ans: a double-type number. 
%--------------------------------------------------------------------------
% 
% Last modified: Jan 29, 2019
% Author: Armin Khamoshi
%-=========================================================================
% Checks to see if the elements of 'list' are unique:
if nargin > 2 
    X(unique(list)) = [];
end

% Find M and allocates memory for the algorithm:
M = length(X);
S = zeros(M, N + 1);

% The summation algorithm:
S(:, 0 + 1) = 1;
S(1, 1 + 1) = X(1);
for i = 2:M
    for j = max(1, i + N - M): min(i, N)
        S(i, j + 1) = S(i - 1, j + 1) + X(i)*S(i - 1, j);
    end
end

% The output:
output = S(M, N+1);
end 
