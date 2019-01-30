function output = AGP_overlap(x, N, list)
%==========================================================================
% This function caluates an element of the n-pair RDM, i.e.
%  <AGP|Pdager_{p1}... Pdager_{pn} P_{q1}...P_{qn}|AGP>. 
%
% Without the ‘list’ input, this function gives <AGP|AGP>.
%
% Inputs:------------------------------------------------------------------
%   x: array of AGP geminals; a Mx1 array.
%   
%   N: the number of pairs in the system. An integer.
%
%   list: a 1D array of 1x2n size. By the definition above, 
%   'list' = [p1,...,pn,q1,...,qn]. Clearly if {p} or {q} have
%   repeated indicies (among themselves), this function gives zero. 
%--------------------------------------------------------------------------
%
% Author: Armin Khamoshi
% Last modified: Dec  2, 2018
%==========================================================================
% Checks and organizes the inputs
if nargin > 2
    l = length(list);
    % 'list' must have even indicies; otherwise it's zero.
    if mod(l,2) 
        output = 0;
        return
    else
        % If the set of {P} or {Pdag} have repetitive indices,
        % give 0.
        l = l/2;
        if l ~= length(unique(list(1:l))) || ...
           l ~= length(unique(list(l+1:end)))
            output = 0;
            return
        end
        k = N - l;
    end
else
    k = N;
    list = [];
end

% To be multiplied later by the overlap:
multiplier = factorial(N)^2*prod(x(list));  

% Prepares the right inputs for ‘sumESF’: 
x(unique(list)) = [];
x = x.^2;

% The 'sumESF' algorithm (calculates the overalp):
M = length(x);
S = zeros(M, k + 1);
S(:, 0 + 1) = 1;
S(1, 1 + 1) = x(1);
for i = 2:M
    for j = max(1, i + k - M): min(i, k)
        S(i, j + 1) = S(i - 1, j + 1) + x(i)*S(i - 1, j);
    end
end

% Final answer:
output = multiplier*S(M, k+1);
end
