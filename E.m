function res = E(x, v, N, R, dynamics)

res = sqrt(B(v, v, N));
% res = B(v, v, N);

% f = @(r) dynamics.cutoff(sqrt(2.*N.*r));
f = @(r) dynamics.a(sqrt(2*N)*r);


lowerbound = sqrt(B(x, x, N));
upperbound = 1000;
temp = integral(f, lowerbound, upperbound);


res = res-  temp;

end

