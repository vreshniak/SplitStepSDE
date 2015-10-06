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
            NPaths  = 1000;
            dt_ref  = 1e-8;
            lev_ini = 16;
            lev_fin = 2^floor(log2(0.1/dt_ref));
            dt_min  = lev_ini  * dt_ref;
            dt_max  = lev_fin * dt_ref;
        case 6
            dt_ref  = 1e-8;
            lev_ini = 16;
            lev_fin = 2^floor(log2(0.1/dt_ref));
            dt_min  = lev_ini * dt_ref;
            dt_max  = lev_fin * dt_ref;
        case 9
            lev_ini = 2^ceil(log2(dt_min/dt_ref));
            lev_fin = 2^floor(log2(dt_max/dt_ref));
            dt_min  = lev_ini * dt_ref;
            dt_max  = lev_fin * dt_ref;
    end

    tspan_ref      = 0:dt_ref:FTime;             % reference time discretization
    tspan_ref(end) = tspan_ref(end-1) + dt_ref;

    M     = size(DiffusionMatrix(1,Y0),2);  % dimension of the noise
    K_ref = length(tspan_ref);              % number of time points
    Nm    = length(method);                 % number of methods to test
    Nlev  = log2(dt_max/dt_min) + 1;        % number of testing levels
       
    err = zeros(Nm,Nlev);
    h   = zeros(1,Nlev);
    
    
    for path = 1:NPaths  
        % reference solution on the fine grid
        Wiener_ref = BrownianMotion(dt_ref,M,K_ref);
        Y_ref = Milstein(DriftVector,DiffusionMatrix,tspan_ref,Y0,Wiener_ref);

        % run test levels
        for p = 1:Nlev
            lev_cur = lev_ini*2^(p-1);
            dt = lev_cur * dt_ref;
            tspan = 0:dt:FTime;
            K = length(tspan);
            Wiener = zeros(M,K);
            for i = 1:K
                j = (i-1)*lev_cur + 1;
                Wiener(:,i) = Wiener_ref(:,j);
            end
    
            % test different methods
            for n = 1:Nm
                Y = method{n}(DriftVector,DiffusionMatrix,tspan,Y0,Wiener);
                err(n, p) = err(n, p) + norm(Y_ref(:,end)-Y(:,end));
            end
            h(p) = dt;
    
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