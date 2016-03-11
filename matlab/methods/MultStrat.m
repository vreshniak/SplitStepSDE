function StratInt = MultStrat(dt,M,dW,p)

%   Purpose
%   =======
%   Find values of the multiple Stratonovich integrals of the form:
%                   /t  /s    j1   j2
%        J(j1,j2) = |   |   dW   dW
%                   /0  /0    
%
%   Method
%   ======
%   Karhunen-Loeve (Fourier) series expansion:
%       Ref - P.Kloeden "Numerical solution of stochastic differential
%                        equation", Chapter 5.8, Chapter 10.3
%
%   Additionally we may utilize the property:
%
%       J(j1,j2) + J(j2,j1) = J(j1)*J(j2)
%           Ref - P.Kloeden "Numerical solution of stochastic differential
%                            equation", (10.3.15)
%
%
%   IN
%   ==
%   1) dt  - integrating time step
%   2) M   - dimension of the white noise
%   3) dW  - increments of driving Wiener processes
%   4) p   - (optional) number of terms in series expansion
%
%   OUT
%   ===
%   StratInt - M-by-M matrix with multiple Ito integrals
%
    
    sqrt2   = sqrt(2);
    sqrt_dt = sqrt(dt);

    % initialize matrix
    StratInt = zeros(M,M);

    if (M > 1)
        % find number of partial sums required to get 1-st order convergence
        % Ref - Kloeden, (10.3.10)   
        if ( p < 0 )
            p = ceil(1/dt);  
        else
            p = min(p,ceil(1/dt));
        end
    
        % additional coefficient
        ro = 0;
        for r = p:-1:1
            ro = ro + 1/(r*r);
        end
        ro = sqrt(1/12 - 0.5*ro/(pi*pi));
    
        mu = sqrt_dt * randn(M,1);
    
        for r = p:-1:1
            nu   = sqrt_dt * randn(M,1);
            zeta = sqrt_dt * randn(M,1);
            for i = 1:M
                StratInt(:,i) = StratInt(:,i) + ( zeta(i)*(sqrt2*dW    + nu   ) - ...
                                                  zeta   *(sqrt2*dW(i) + nu(i)) ) / r;
            end
        end
        
        for i = 1:M
            StratInt(:,i) = StratInt(:,i)/pi + 2*ro*(mu(i)*dW - mu*dW(i));
        end
    end
    
    for i = 1:M
        StratInt(:,i) = StratInt(:,i) + dW(i)*dW;
    end
    StratInt = 0.5 * StratInt;

end



























