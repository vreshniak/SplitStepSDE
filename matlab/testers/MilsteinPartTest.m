function MilsteinPartTest()

%   Purpose
%   =======
%   Test computational complexity and relative accuracy of different
%   variants of approxiamtion of the Milstein part of the numerical 
%   method
%


    M     = 100;            % stochastic dimension
    N     = 100;            % number of equations
    Ito   = randn(M,M);     % Ito matrix
    G     = randn(N,M);     % diffusion matrix
    dg_dx = randn(N,M,N);   % DiffusionMatrixGradient
    dg_dx_T = zeros(M,N,N); % transposed DiffusionMatrixGradient
    for i = 1:N
        dg_dx_T(:,:,i) = dg_dx(:,:,i)';
    end

    for i = 1:10
        tic; MilPart_1  = MilsteinPart1();  t1 = toc;
        tic; MilPart_2a = MilsteinPart2a(); t2 = toc;
        tic; MilPart_2b = MilsteinPart2b(); t3 = toc;
        tic; MilPart_3  = MilsteinPart3();  t4 = toc;
    end
    
    err1 = norm(MilPart_2a - MilPart_1);
    err2 = norm(MilPart_2b - MilPart_1);
    err3 = norm(MilPart_3  - MilPart_1);
    
    fprintf('Algorithm 1 : %e seconds\n',t1);
    fprintf('Algorithm 2a: %e seconds\n',t2);
    fprintf('Algorithm 2b: %e seconds\n',t3);
    fprintf('Algorithm 3 : %e seconds\n',t4);
    fprintf('\n');
    fprintf('Error between algorithms 1 and 2: %e\n',err1);
    fprintf('Error between algorithms 1 and 3: %e\n',err2);
    fprintf('Error between algorithms 1 and 4: %e\n',err3);

    
   
	function result = MilsteinPart1()
        result = zeros(N,1);
        dg_dx_1 = permute(dg_dx,[1 3 2]);
        B = G * Ito;
        for j = 1:M
            result = result + dg_dx_1(:,:,j)*B(:,j);
        end
    end

    function result = MilsteinPart2a()
        % Ref. Reshniak et al. 
        %   "Split-step Milstein methods for multi-channel stiff stochastic 
        %    differential systems". Applied Numerical Mathematics
        %  formula (A.3)
        result = zeros(N,1);
        B = (G * Ito)';
        for j = 1:N
            result = result + dg_dx(:,:,j)*B(:,j);
        end
    end

    function result = MilsteinPart2b()
        % Ref. Reshniak et al. 
        %   "Split-step Milstein methods for multi-channel stiff stochastic 
        %    differential systems". Applied Numerical Mathematics
        %  formula (A.3)
        result = zeros(N,1);
        B = G * Ito;
        for j = 1:N
            result = result + (B(j,:)*dg_dx_T(:,:,j))';
        end
    end

    function result = MilsteinPart3()
        % classical triple sum
        result = zeros(N,1);
        for j1 = 1:M
            for j2 = 1:M
                for k = 1:N
                    result = result + dg_dx(:,j2,k)*G(k,j1)*Ito(j1,j2);
                end
            end
        end
    end

end
