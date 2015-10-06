function [G,N,M] = DiffusionMatrix_apnum2015_ex2(t,X,transp)

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


    N = 2;
    M = 2;

    if ( transp > 0 ) % transposed matrix
        G = [ X(1) + X(2),   X(1) + X(2); ...
              X(1) - X(2),  -X(1) + X(2) ];
    else
        G = [ X(1) + X(2),   X(1) - X(2); ...
              X(1) + X(2),  -X(1) + X(2) ];
    end

end



















