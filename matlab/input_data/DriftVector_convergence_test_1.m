function [F,N] = DriftVector_convergence_test_1(t,X)

%   Purpose
%   =======
%   Calculate drift vector
%
%   IN
%   ==
%   1) t - time
%   2) Y - N-by-1 vector of solution at time t
%
%   OUT
%   ===
%   F - N-by-1 drift vector
%   N - dimension of the system


    N = 1;

    a = 1.5;
    
	F = a * X(1);
   
end















