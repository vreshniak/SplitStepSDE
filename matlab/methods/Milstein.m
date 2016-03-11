function [Y,Wiener] = Milstein(DriftVector,DiffusionMatrix,T,Y0,Wiener,p)

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
%                               __M                    __M    __M    __N
%                               \                      \      \      \      dg(i,j2)
%   yi[k+1] = yi[k] + h*fi[k] + /     g(i,j) I(j)   +  /      /      /     --------- g(k,j1) I(j1,j2)
%                               --j=1                  --j1=1 --j2=1 --k=1    dyk
%                                               
%
%   IN
%   ==
%   1) DriftVector     - function handle that evaluates drift vector
%   2) DiffusionMatrix - function handle that evaluates matrix of diffusion coefficients
%   3) T  - K-dimensional vector of time points
%   4) Y0 - N-dimensional column vector with initial data
%   5) Wiener - (optional) M-by-K array of driving Wiener processes (same as in the output below)
%   4) p      - (optional) number of terms in series expansion of multiple stochastic integrals
%
%
%   OUT
%   ===
%   Y      - N-by-K solution array. Each row in Y is the solution of the corresponding equation
%   Wiener - M-by-K array of the driving Wiener processes. 


    % number of equations and dimension of the noise
    [~,N,M] = DiffusionMatrix(1,Y0);
    
    % number of points in time discretization
    K = length(T);

    % step size
    dt = T(2) - T(1);

    % initialize array with gradient of the diffusion matrix
    dg_dx = zeros(N,M,N);
    
    % initialize solution array
    Y = zeros(N,K);
    Y(:,1) = Y0(:);
    
    % generate array of driving Wiener processes
    switch nargin 
        case 4
            Wiener = BrownianMotion(dt,M,K);
            p = 50;
        case 5
            p = 50;
    end
    
    % loop in time
    for i = 2:K
        % generate vector of noise increments
        dW = Wiener(:,i) - Wiener(:,i-1);

        F = DriftVector(T(i-1),Y(:,i-1));
        G = DiffusionMatrix(T(i-1),Y(:,i-1));
        
        Ito = MultIto(dt,M,dW,p);
                
        % update solution
        Y(:,i) = Y(:,i-1) + F*dt + G*dW;  addMilsteinPart();
    end
      
    
    
    function addMilsteinPart()
        % Ref. Reshniak et al. 
        %   "Split-step Milstein methods for multi-channel stiff stochastic 
        %    differential systems". Applied Numerical Mathematics
        %  formula (A.3)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find gradient of the diffusion matrix (see DiffusionMatrixGradient)
        dx = 1e-6;
        Y1 = Y(:,i-1);
        for k = 1:N
            Y1(k) = Y(k,i-1) + dx;
            dg_dx(:,:,k) = DiffusionMatrix(T(i-1),Y1) - G; 
            Y1(k) = Y(k,i-1);
        end
        dg_dx = dg_dx / dx; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        B = (G*Ito)';
        for j = 1:N
            Y(:,i) = Y(:,i) + dg_dx(:,:,j)*B(:,j);
        end
    end

end


























