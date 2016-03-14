function [F,G,X0] = apnum2015_ex3()
%   Purpose
%   =======
%
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


	% reaction rate constants
	c = [ 1.0d3;  ...
          1.0d3;  ...
          1.0d-5; ...
          1.0d1;  ...
          1.0d0;  ...
          1.0d6 ];
      
    % stochiometric coefficients
    nu = [ -1  1 -1  1  1 -1; ...
           -1  1  1 -1 -1  1; ...
            1 -1 -1  1 -1  1];

    F  = @DriftVector;
    G  = @DiffusionMatrix;
	X0 = [1e3; 1e3; 1e6];

    
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

        N = 3;
    
        F = zeros(N,1);
    
    	% propensities
    	alpha = [ c(1) * X(1) * X(2); ...
                  c(2) * X(3);        ...
                  c(3) * X(1) * X(3); ...
                  c(4) * X(2);        ...
                  c(5) * X(2) * X(3); ...
                  c(6) * X(1) ];
    
        for i = 1:N
            for j = 1:6
                F(i) = F(i) + nu(i,j) * alpha(j);
            end
        end 
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
    
        N = 3;
        M = 6;
    
        G = zeros(N,M);
    
    	% propensities
    	alpha = [ c(1) * X(1) * X(2); ...
                  c(2) * X(3);        ...
                  c(3) * X(1) * X(3); ...
                  c(4) * X(2);        ...
                  c(5) * X(2) * X(3); ...
                  c(6) * X(1) ];
    
        for i = 1:N
            for j = 1:M
                G(i,j) = nu(i,j) * sqrt(alpha(j));
            end
        end
    end
    
end




