function [Y,Wiener] = Milstein(DriftVector,DiffusionMatrix,T,Y0,Wiener)

%   Purpose
%   =======
%   Find solution of the system of Ito stochastic equations with 
%   multi-channel non-commutative noise:
%
%      / Y1 \   / f1 \        / g11  g12 ... g1m \   / dW1 \
%      | Y2 |   | f2 |        | g12  g22 ... g2m |   | dW2 |
%      | .  | = | .  | * dt + |  .      .     .  | * |  .  |
%      | .  |   | .  |        |  .         .  .  |   |  .  |
%      \ Yn /   \ fn /        \ gn1  gn2 ... gnm /   \ dWm /
%
%       Yi(t0) = Yi0, i = 1..n
%
%
%   Method
%   ======
%   Milstein method on uniform time grid:
%                               __M             
%                               \               
%   yi[k+1] = yi[k] + h*fi[k] + /     g(i,j) I(j)
%                               --j=1           
%                                               
%                      __M    __M    __N
%                      \      \      \      dg(i,j2)
%                   +  /      /      /     --------- g(k,j1) I(j1,j2)
%                      --j1=1 --j2=1 --k=1    dyk
%
%
%   IN
%   ==
%   1) DriftVector     - function handle that evaluates drift vector
%   2) DiffusionMatrix - function handle that evaluates matrix of 
%                        diffusion coefficients
%   3) T  - vector of time points
%   4) Y0 - N-dimensional column vector with initial data
%   5) Wiener - optional array of driving Wiener processes 
%                 (same as in the output below)
%
%
%   OUT
%   ===
%   Y - N-by-K solution array. Each row in Y is the solution of the
%       corresponding equation
%   Wiener - M-by-K-dimensional array of the driving Wiener processes. 


    % number of equations and dimension of the noise
    [G,N,M] = DiffusionMatrix(1,Y0);
    
    % number of points in time discretization
    K = length(T);

    % step size
    dt = T(2) - T(1);

    % initialize solution array
    Y = zeros(N,K);
    Y(:,1) = Y0(:);
    
    % generate array of driving Wiener processes
    if ( nargin == 4 )
        Wiener = BrownianMotion(dt,M,K);
    end
    
    % loop in time
    for i = 2:K
        % generate vector of noise increments
        dW = Wiener(:,i) - Wiener(:,i-1);

        F = DriftVector(T(i-1),Y(:,i-1));
        G = DiffusionMatrix(T(i-1),Y(:,i-1));
        
        Ito = MultIto(dt,M,dW);
                
        % update solution
        Y(:,i) = Y(:,i-1) + F*dt + G*dW + MilsteinPart(T(i-1),Y(:,i-1));
    end
      
    
    function result = MilsteinPart(t,X)
        % Ref. Reshniak et al. 
        %   "Split-step Milstein methods for multi-channel stiff stochastic 
        %    differential systems". Applied Numerical Mathematics
        %  formula (A.3)
        result = zeros(N,1);
        dg_dx = DiffusionMatrixGradient(t,X,N,M);
        B = (G*Ito)';
        for j = 1:N
            result = result + dg_dx(:,:,j)*B(:,j);
        end
    end

end


























