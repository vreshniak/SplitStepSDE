function [h,err,order] = StrongConvergenceTestWithExactSol_1(method,FTime,NPaths,dt_min,dt_max)

%   Purpose
%   =======
%   Test strong convergence of the numerical method using system wich has exact solution. 
%
%   Method
%   ======
%   Strong error is obtained using the exact solution of the system
%
%              dX(t) = a*X(t)dt + b*X(t)dW(t)
%               X(0) = 1
%
%   IN
%   ==
%   1) method - cell array of function handles that calls corresponding numerical methods
%   2) FTime  - final time of the simulation
%   3) Y0     - n-dimensional column vector with initial data
%   4) NPaths - (optional) number of samples for Monte Carlo estimator
%   5) dt_ref - (optional) reference time discretization step
%   6) dt_min - (optional) finest time discretization step
%   7) dt_max - (optional) coarsest time discretization step
%
%
%   OUT
%   ===
%   h     - vector with time steps
%   err   - vector with corresponding strong errors
%   order - estimated orders of convergence


    switch nargin 
        case 2
            NPaths = 100;
            dt_min = 2^(-15);
            dt_max = 2^(-3);
        case 3
            dt_min = 2^(-15);
            dt_max = 2^(-3);
        case 5
%             dt_min = 2^ceil(log2(dt_min));
            dt_max = 2^floor(log2(dt_max/dt_min)) * dt_min;
    end


    fprintf('dt_min = %6.4e\n', dt_min);
    fprintf('dt_max = %6.4e\n', dt_max);
    
    
    M     = 1;                              % dimension of the noise
    Nm    = length(method);                 % number of methods to test
    Nlev  = log2(dt_max/dt_min) + 1;        % number of testing levels
       
    % allocate various arrays
    err      = zeros(Nm,Nlev);              % strong errors
    h        = zeros(1,Nlev);               % time steps 
    grid_cur = zeros(1,Nlev);               % levels of discretization
    tspan    = cell(1,Nlev);                % time points
    K        = zeros(1,Nlev);               % number of time points
    Wiener   = cell(1,Nlev);                % Brownian motions
	for p = 1:Nlev
        grid_cur(p) = 2^(p-1);
        h(p)        = grid_cur(p) * dt_min;
        tspan{p}    = 0:h(p):FTime;
        K(p)        = length(tspan{p});
        Wiener{p}   = zeros(M,K(p));
	end
    
    % parameters of the stochastic system
    a  = 1;
    b  = 1;
    Y0 = 1;
    
    for path = 1:NPaths  
        % reference solution on the fine grid
        Wiener_ref = BrownianMotion(dt_min,M,K(1));
        Y_ref = Y0*exp( (a-0.5*b^2)*FTime + b*Wiener_ref(end) );

        % run test levels
        for p = 1:Nlev
            % extract "lower dimensional" Brownian motion
            for i = 1:K(p)
                Wiener{p}(:,i) = Wiener_ref(:,(i-1)*grid_cur(p)+1);
            end
    
            % test different methods
            for n = 1:Nm
                Y = method{n}(@DriftVector,@DiffusionMatrix,tspan{p},Y0,Wiener{p});
                err(n, p) = err(n, p) + norm(Y_ref-Y(:,end));
            end
    
            fprintf('Path: %i of %i    Level: %i\n',path,NPaths,p);
        end
    end
   
    err = err ./ NPaths;
    
    order = zeros(Nm, 1);
    for n = 1:Nm
        pp = polyfit( log(h), log(err(n,:)), 1 );
        order(n) = pp(1);
    end
    
    
    function [F,N] = DriftVector(t,X)
        N = 1;
    	F = a * X(1);
    end

    function [G,N,M] = DiffusionMatrix(t,X)
        N = 1;
        M = 1;
        G = b * X(1);
    end
    

end