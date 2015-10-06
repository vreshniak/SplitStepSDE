function [F,N] = DriftVector_apnum2015_ex2(t,X)

%   Purpose
%   =======
%   Calculate drift vector
%
%   Reference:
%   Reshniak, V., Khaliq, A. Q. M., Voss, D. A., & Zhang, G. (2015). 
%   "Split-step Milstein methods for multi-channel stiff stochastic 
%   differential systems". Applied Numerical Mathematics
%
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


    N = 2;

	F = [ X(2); ...
         -X(1) ];
   
end















