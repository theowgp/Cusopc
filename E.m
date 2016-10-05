function res = E(x, v, N, R, dynamics)

res = sqrt(B(v, v, N));
% res = B(v, v, N);

% f = @(r) dynamics.cutoff(sqrt(2.*N.*r));
f = @(r) tempf(sqrt(2*N*r), R);


lowerbound = sqrt(B(x, x, N));
upperbound = R^2/(2*N);
temp = integral(f, lowerbound, upperbound);


res = res-  temp;

end

