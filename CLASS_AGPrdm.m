% This class facilitates most calculations concerning AGP RDM.
% 
% To see the necessary inputs for creating this class, see the constructor
% function.
%
% To see the analytic expressions and conventions used in this function,
% please refer to the manual.
%
%--------------------------------------------------------------------------
% public functions: (can be accessed by instantiating the class)
%
%   print_AGP_RDM: creates .mat files containing the RDM's associated with 
%                  each value of G. ATTENTION: can only use this function 
%                  when you are in the directory that contains the eta’s.
%-------------------------------------------------------------------------- 
% static functions: (can be accessed without instantiating the class)
%
%   nRDM: This function caluates the n-RDM over AGP.   
%
%   test_myRDM: Sees if an RDM is correctly calculated by random sampling.               
%
%   UniqueIndices_AGPrdm: Creates an array of unique indices associated with 
%                         the elements of n-body RDM over AGP.
%
%   sumESF: performs the sum of elementary symmetric function.
%   UniqueIndices_AGPrdm: Creates the indices needed in computing the RDM.
% ------------------------------------------------------------------------   
% 
% Last modified: Jan 29, 2019
% Author: Armin Khamoshi
%==========================================================================
classdef CLASS_AGPrdm < handle
    properties (SetAccess = 'protected')   
        numLevels            % Number of levels
        numPairs             % Number of pairs
        
        % Values of G for which we have data:
        G = [0.1 0.3 0.500:0.1:2.00 2.50 3.00]; 
    end
    properties (Hidden)
        % To contain the indices for 2, 4, 6 RDM's:
        uniqueIndices = cell(3,1);  
        
        % Will be used to impose symmetries of the creation and 
        % annihilation operators on the RDM:
        symm_4RDM = [ 2 1 3 4; 
                      1 2 4 3
                      3 4 1 2];
        symm_6RDM = [ 4 5 6 1 2 3;
                      2 1 3 4 5 6;
                      3 2 1 4 5 6;
                      1 3 2 4 5 6;
                      1 2 3 5 4 6;
                      1 2 3 6 5 4;
                      1 2 3 4 6 5];
    end
    methods
        function obj = CLASS_AGPrdm(numPairs, numLevels, n)
            %==============================================================
            % The constructor of the 'CLASS_AGPrdm' class:
            %
            % Inputs:------------------------------------------------------ 
            %   numLevels: the number of orbitals. An inetger.   
            %   numPairs: the number of pairs in the system. An integer.
            %   n: (optiona) the dimensionality of the RDM. (e.g. 2, 4, 6).
            % -------------------------------------------------------------
            %
            % Last modified: Jan 21, 2019
            % Author: Armin Khamoshi
            %==============================================================
            obj.numLevels = numLevels;
            obj.numPairs = numPairs; 
            if nargin > 2
                obj.print_AGP_RDM(n);
            end
        end
        function print_AGP_RDM(obj, n)
            %==============================================================
            % Reads the files containing eta’s, calculates the RDM, and
            % prints .mat files containing the RDM associated with each.
            %
            % The names of files containing the eta’s must be stores as
            % follows: eta_[number of levels]L[number of pairs]P_[the value
            % of G].mat. (Without the brackets, obviously.)
            %
            % ATTENTION: can only use this function then you are in the  
            % directory that contains the eta’s.
            %
            % The list of values of G is stored in obj.G (see the
            % properties.)
            %
            % Inputs:------------------------------------------------------
            %   n: the rank of the RDM
            % Outputs:-----------------------------------------------------
            %   .mat files containing the RDM associated with each value of
            %   G.
            % -------------------------------------------------------------
            %
            % Last modified: Jan 28, 2019
            % Author: Armin Khamoshi
            %==============================================================
            M = obj.numLevels;
            N = obj.numPairs;
            
            % Checks to see if the the indices have already been calculated
            if isempty(obj.uniqueIndices{n/2})
                obj.uniqueIndices{n/2} = obj.UniqueIndices_AGPrdm(n, M);
            end
            
            % Allocates memory for the RDM:
            RDM = zeros(M*ones(1, n));
            
            for i = obj.G
                %imports eta’s the current directory
                eta = importdata(['eta_', num2str(M), 'L', num2str(N), 'P_', num2str(i), '.mat']);

                % Populates the 'necessary' elements only: 
                for j = 1:length(obj.uniqueIndices{n/2})
                    list = obj.uniqueIndices{n/2}(j,:);
                    arg = num2cell(list, 1);
                    RDM(arg{:}) = obj.nRDM(eta, N, list);
                end

                % Uses the symmetries to find the remaining matrix elements 
                switch n
                    case 2
                        RDM = max(RDM, RDM');
                    case 4
                        for j = 1:size(obj.symm_4RDM,1)
                            RDM = max(RDM, permute(RDM, obj.symm_4RDM(j,:)));
                        end
                    case 6
                        for j = 1:size(obj.symm_6RDM,1)
                            RDM = max(RDM, permute(RDM, obj.symm_6RDM(j,:)));
                        end
                end

                % Exports the RDM by saving it as a .mat file.
                numbers = {'two', 'four', 'six', 'eight'}; 
                filename = [numbers{n/2}, 'RDM_', num2str(M), 'L', ...
                            num2str(N), 'P_', num2str(i), '.mat'];
                save(filename, 'RDM');
                disp(['Printing the RDM for G = ', num2str(i), ' is finished.']);
            end            
        end    
    end
    methods (Static)
        function output = sumESF(X, N, list)
            %==============================================================
            % Calculates the Nth Elementary Symmetric Function (ESF) 
            %  associated with a vector of M numbers X = ( x_1,...,x_M )
            %  such that the indecies in array ‘list' are eliminated. 
            % 
            % The algorithm for carrying the sum is copied from:
            % “Accurate validated and fast evaluation of elementary 
            % symmetric functions and its application,” Hao Jiang et al.
            %
            % Inputs:------------------------------------------------------
            %   X: is a 1xM dimentional array as explained above.
            %   N: is the number of products in each summand. Must be an 
            %   integer. list: contains the indecies that are to be
            %   eliminates..
            %--------------------------------------------------------------
            %
            % Author: Armin Khamoshi
            % Last modified: Nov  21, 2018
            %==============================================================
            % Checks the inputs
            if nargin > 2 
                X(unique(list)) = [];
            end
            
            % The 'sumESF' algorithm:
            M = length(X);
            S = zeros(M, N + 1);
            S(:, 0 + 1) = 1;
            S(1, 1 + 1) = X(1);
            for i = 2:M
                for j = max(1, i + N - M): min(i, N)
                    S(i, j + 1) = S(i - 1, j + 1) + X(i)*S(i - 1, j);
                end
            end
                output = S(M, N+1);
        end
        function output = nRDM(x, N, list)
            %==============================================================
            % This function caluates the n-RDM, i.e.
            %  <AGP|Pdager_{p1}... Pdager_{pn} P_{q1}...P_{qn}|AGP>. 
            %
            % Without any inputs, this function gives <AGP|AGP>.
            %
            % Inputs:------------------------------------------------------
            %   x: array of AGP geminals; a Mx1 array.
            %
            %   N: the number of pairs in the system. An integer.
            %
            %   list: a 1D array of 1x2n size. By the definition above, 
            %   'list' = [p1,...,pn,q1,...,qn]. Clearlif ~isempty
            %   (uniqueIndices_2RDM)y if {p} or {q} have repeated indicies,
            %   this function gives zero. 
            % -------------------------------------------------------------
            %
            % Author: Armin Khamoshi
            % Last modified: Jan  29, 2018
            %==============================================================
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
        function output = UniqueIndices_AGPrdm(n, M)
            %==============================================================
            % Produces unique indices associated with the n-body RDM over AGP.
            % Explicitly, let g_pqrs = <AGP| Pdag_p Pdag_q P_r P_s |AGP>. Then
            % this function gives an array containing all permutations of
            % g_pqrs such that
            % p < q, and r < s.
            %
            % Inputs:------------------------------------------------------
            %   n = the dimentionality of RDM; in the example abFortranPrinterove, n = 4.
            %   Must be an even integer.
            %      
            %   M = the number of orbitals. An integer.
            % Output:------------------------------------------------------
            %   output = unique indices associated with the n-body RDM over
            %   AGP. Has dimentions 1/2*M_choose_n*(M_choose_n + 1) by n.  
            %--------------------------------------------------------------
            %
            % Author: Armin Khamoshi
            % Last modified: Jan  17, 2019
            %==============================================================
            n = n/2;    % For convenience:       

            % Creates an array containing all unique indices of the string  
            % of creation/annihilation operators:
            Creation_indices = nchoosek(1:M, n);
            Annihilation_indices = Creation_indices;

            % Allocates memory for all unique indices:
            M_choose_n = nchoosek(M, n);
            Allpermutaions = zeros(1/2*M_choose_n*(M_choose_n + 1) , 2*n);

            % Attaches the indices of the creation operators to those of the 
            % annihilation operator:
            k = 0;
            for  i = 1:M_choose_n 
                for j = i:M_choose_n
                    k = k + 1;
                    Allpermutaions(k,:) = [Creation_indices(i,:), ...
                        Annihilation_indices(j,:)];
                end
            end

            % Produces the output
            output = Allpermutaions;
        end
        function test_myRDM(size, N, eta, RDM)
            % =============================================================
            % This function checks if the RDM has been calculated correctly. 
            % It does so by randomly sampling certain indices of the RDM and 
            % calculating the corresponding element again on the fly to see 
            % of they match. 
            %
            % Inputs:------------------------------------------------------
            %   size: The size of the sample. A positive integer.
            %   N: the number of levels
            %   eta: the array of geminals. Mx1 array.
            %   RDM: the RDM matrix.
            % -------------------------------------------------------------
            %
            % Author: Armin Khamoshi
            % Last updated: Jan 28, 2019
            % =============================================================
            % Preliminary stuff:
            M = length(eta);    % number of levels
            n = ndims(RDM);     % dimensionality of the RDM
            erFound = 0;        % To signal of an error has been found
            
            for i = 1:size
                % Generates a random index toeck the RDM with:
                sample = randperm(M, n);
                
                % Changes the data-type of 'sample' for later use:
                arg = num2cell(sample, 1);
                
                % Print where the inconstancy has been found (if at all):
                if ~(abs(CLASS_AGPrdm.nRDM(eta, N, sample) - RDM(arg{:}) )...
                        < 1E-15) 
                    disp(['Found an inconstancy at: ', num2str(sample)]);
                    erFound = 1;
                end
            end
            
            % If no inconsistencies found...
            if ~erFound 
                disp('No inconsistencies found in this random ensemble.')
            end
        end
    end 
end

