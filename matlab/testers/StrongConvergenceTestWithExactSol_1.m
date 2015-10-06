function [h,err,order] = StrongConvergenceTestWithExactSol_1(method,FTime,NPaths,dt_ref,dt_min,dt_max)

%   Purpose
%   =======
%   Test strong convergence of the numerical method using system wich 
%   has exact solution. 
%   Noise channels are identical -> the simplest case of commutative noise
%
%   Reference:
%   Reshniak, V., Khaliq, A. Q. M., Voss, D. A., & Zhang, G. (2015). 
%   "Split-step Milstein methods for multi-channel stiff stochastic 
%   differential systems". Applied Numerical Mathematics
%
%
%   Method
%   ======
%   Strong error is obtained using the exact solution
%
%   IN
%   ==
%   1) method          - function handle that calls the corresponding
%                        numerical method
%   2) FTime - final time of the simulation
%   3) NPaths - (optional) number of samples for Monte Carlo estimator
%   4) dt_ref - (optional) reference time discretization step
%   5) dt_min - (optional) finest time discretization step
%   6) dt_max - (optional) coarsest time discretization step
%
%
%   OUT
%   ===
%   h     - vector with time steps
%   err   - vector with corresponding strong errors
%   order - estimated order of convergence

    switch nargin 
        case 2
            NPaths = 1000;
            dt_ref = 1e-8;
            dt_min = 1e-5;
            dt_max = 1e-1;
        case 3
            dt_ref = 1e-8;
            dt_min = 1e-5;
            dt_max = 1e-1;
    end
    
	lev_ini = 2^ceil(log2(dt_min/dt_ref));
    lev_fin = 2^floor(log2(dt_max/dt_ref));
    dt_min  = lev_ini * dt_ref;
    dt_max  = lev_fin * dt_ref;
    
    tspan_ref = 0:dt_ref:FTime;             % reference time discretization
    tspan_ref(end) = tspan_ref(end-1) + dt_ref;

    N     = 5;                              % number of equations
    M     = 20;                             % number of noise channels
    K_ref = length(tspan_ref);              % number of time points
    Nm    = length(method);                 % number of methods to test
    Nlev  = log2(dt_max/dt_min) + 1;        % number of testing levels

    A = ones(N) / 20;   A(logical(eye(N))) = -1.5;
    B = ones(N) / 100;  B(logical(eye(N))) =  0.2;
    
    DriftVector     = @(t,X)A*X;
    DiffusionMatrix = @(t,X)repmat(B*X,1,M);
    
    Y0 = ones(N,1);
    
    err = zeros(Nm,Nlev);
    h   = zeros(1,Nlev);

    for path = 1:NPaths  
        % reference solution on the fine grid
        Wiener_ref = BrownianMotion(dt_ref,M,K_ref);
        Y_ref = zeros(N,K_ref);
        for i = 1:K_ref
            Y_ref(:,i) = expm( ( A - 0.5*M*B*B )*tspan_ref(i) + B*sum(Wiener_ref(:,i)) ) * Y0;
        end

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
    
            fprintf('Path: %i of %i    Level: %i of %i\n',path,NPaths,p,Nlev);
        end

    end

    err = err ./ NPaths;

    order = zeros(Nm, 1);
    for n = 1:Nm
        p = polyfit( log(h), log(err(n,:)), 1 );
        order(n) = p(1);
    end

end