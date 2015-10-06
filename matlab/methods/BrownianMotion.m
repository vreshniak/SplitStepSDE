function B = BrownianMotion(dt,M,K)

%   Purpose
%   =======
%   Generate M independent Wiener processes
%
%   IN
%   ==
%   1) dt  - integrating time step
%   2) M   - dimension of the white noise
%   3) K   - number of time points
%
%   OUT
%   ===
%   B - M-by-K matrix with Wiener processes

    sqrt_dt = sqrt(dt);

	B = zeros(M,K);
    for i = 2:K
        B(:,i) = B(:,i-1) + sqrt_dt * randn(M,1);
    end

end