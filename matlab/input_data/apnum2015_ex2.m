function [F,G,X0] = apnum2015_ex2()
%   Purpose
%   =======
%   Reference:
%   Reshniak, V., Khaliq, A. Q. M., Voss, D. A., & Zhang, G. (2015). 
%   "Split-step Milstein methods for multi-channel stiff stochastic 
%   differential systems". Applied Numerical Mathematics
%
%   OUT
%   ===
%   F  - function handle to drift vector
%   G  - function handle to diffusion matrix
%   Y0 - initial data

    F  = @DriftVector;
    G  = @DiffusionMatrix;
    X0 = [1; 0];
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [F,N] = DriftVector(t,X)
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
    
        N = 2;
    
    	F = [ X(2); ...
             -X(1) ];   
    end


    function [G,N,M] = DiffusionMatrix(t,X)
    %   Purpose
    %   =======
    %   Calculate diffusion matrix
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
    
        G = [ X(1) + X(2),   X(1) - X(2); ...
              X(1) + X(2),  -X(1) + X(2) ];
    end

end


