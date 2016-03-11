function [h,err,order] = StrongConvergenceTest(method,DriftVector,DiffusionMatrix,FTime,Y0,NPaths,dt_ref,dt_min,dt_max)

%   Purpose
%   =======
%   Test strong convergence of the numerical method 
%
%   Method
%   ======
%   Strong error is obtained using the reference solution calculated
%   with the Milstein method on the fine reference grid
%
%   IN
%   ==
%   1) method          - cell array of function handles that calls 
%                        corresponding numerical methods
%   2) DriftVector     - function handle that evaluates drift vector
%   3) DiffusionMatrix - function handle that evaluates matrix of 
%                        diffusion coefficients
%   4) FTime  - final time of the simulation
%   5) Y0     - n-dimensional column vector with initial data
%   6) NPaths - (optional) number of samples for Monte Carlo estimator
%   7) dt_ref - (optional) reference time discretization step
%   8) dt_min - (optional) finest time discretization step
%   9) dt_max - (optional) coarsest time discretization step
%
%
%   OUT
%   ===
%   h     - vector with time steps
%   err   - vector with corresponding strong errors
%   order - estimated orders of convergence


    switch nargin 
        case 5
            NPaths   = 1000;
            dt_ref   = 1e-8;
            grid_min = 16;
            grid_max = 2^floor(log2(0.1/dt_ref));
            dt_min   = grid_min * dt_ref;
            dt_max   = grid_max * dt_ref;
        case 6
            dt_ref   = 1e-8;
            grid_min = 16;
            grid_max = 2^floor(log2(0.1/dt_ref));
            dt_min   = grid_min * dt_ref;
            dt_max   = grid_max * dt_ref;
        case 9
            grid_min = 2^ceil(log2(dt_min/dt_ref));
            grid_max = 2^floor(log2(dt_max/dt_ref));
            dt_min   = grid_min * dt_ref;
            dt_max   = grid_max * dt_ref;
    end

    % reference time discretization
    tspan_ref      = 0:dt_ref:FTime;        
    tspan_ref(end) = tspan_ref(end-1) + dt_ref;

    M     = size(DiffusionMatrix(1,Y0),2);  % dimension of the noise
    K_ref = length(tspan_ref);              % number of time points
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
        grid_cur(p) = grid_min*2^(p-1);
        h(p)        = grid_cur(p) * dt_ref;
        tspan{p}    = 0:h(p):tspan_ref(end);
        K(p)        = length(tspan{p});
        Wiener{p}   = zeros(M,K(p));
	end
    
    for path = 1:NPaths  
        % reference solution on the fine grid
        Wiener_ref = BrownianMotion(dt_ref,M,K_ref);
        Y_ref = Milstein(DriftVector,DiffusionMatrix,tspan_ref,Y0,Wiener_ref);
%         Y_ref = EulerMaruyama(DriftVector,DiffusionMatrix,tspan_ref,Y0,Wiener_ref);

        % run test levels
        for p = 1:Nlev
            % extract "lower dimensional" Brownian motion
            for i = 1:K(p)
                Wiener{p}(:,i) = Wiener_ref(:,(i-1)*grid_cur(p)+1);
            end
    
            % test different methods
            for n = 1:Nm
                Y = method{n}(DriftVector,DiffusionMatrix,tspan{p},Y0,Wiener{p});
                err(n, p) = err(n, p) + norm(Y_ref(:,end)-Y(:,end));
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

end