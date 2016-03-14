function [F,G,X0] = convergence_test_1()
%   Purpose
%   =======
%
%   OUT
%   ===
%   F  - function handle to drift vector
%   G  - function handle to diffusion matrix
%   Y0 - initial data

    F  = @DriftVector;
    G  = @DiffusionMatrix;
	X0 = 1.0;

    
    function [F,N] = DriftVector(t,X)
    %   Purpose
    %   =======
    %   Calculate drift vector
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

        N = 1;
    
        a = 1.5;
    
    	F = a * X(1);
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
    
        N = 1;
        M = 1;
    
        b = 1;
    
        G = b * X(1);
	end
    
end




