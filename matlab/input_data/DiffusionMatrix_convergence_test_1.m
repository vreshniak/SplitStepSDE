function [G,N,M] = DiffusionMatrix_convergence_test_1(t,X)

%   Purpose
%   =======
%   Calculate diffusion matrix
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
%   2) Y - N-dimensional vector of solution at time t
%
%   OUT
%   ===
%   G - N-by-M diffusion matrix
%   N - dimension of the system
%   M - dimension of the noise


    N = 1;
    M = 1;
    
    b = 1;

    G = b * X(1);

end



















