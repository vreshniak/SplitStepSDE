    function dg_dx = DiffusionMatrixGradient(t,X,DiffusionMatrix,N,M)
    
    %   Purpose
    %   =======
    %   Find gradient of the diffusion matrix:
    %                                 /-----------------------------------\
    %                                 | dG11/dXn  dG12/dXn  ...  dG1m/dXn |
    %                      /-----------------------------------\      .   |
    %                      | dG11/dX2  dG12/dX2  ...  dG1m/dX2 |      .   |
    %           /-----------------------------------\ dG2m/dX2 | dGnm/dXn |
    %           | dG11/dX1  dG12/dX1  ...  dG1m/dX1 |      .   | ---------/
    %           | dG21/dX1  dG22/dX1  ...  dG2m/dX1 |      .   |     .
    %   dg/dx = |   .                 .         .   | dGnm/dX2 |  .
    %           |   .                    .      .   | ---------/
    %           | dGn1/dX1  dGn2/dX1  ...  dGnm/dX1 |
    %           \-----------------------------------/
    %
    %   IN
    %   ==
    %   1) t - time
    %   2) X - N-by-1 vector of solution at time t
    %
    %   OUT
    %   ===
    %   dg_dx - N-by-M-by-N array. Each n-th N-by-M slice of the array
    %           corresponds to the derivative of the diffusion matrix with
    %           respect to Xn
    
        dx = 1e-7;
        dg_dx = zeros(N,M,N);
        G_ref = DiffusionMatrix(t,X);
        X1 = X;
        for k = 1:N
            X1(k) = X(k) + dx;
            dg_dx(:,:,k) = DiffusionMatrix(t,X1) - G_ref; 
            X1(k) = X(k);
        end
        dg_dx = dg_dx / dx; 
        
    end