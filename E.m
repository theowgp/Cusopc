% function res = E(x, v, N, R, dynamics)
% 
% res = sqrt(B(v, v, N));
% % res = B(v, v, N);
% 
% % f = @(r) dynamics.cutoff(sqrt(2.*N.*r));
% f = @(r) dynamics.a(sqrt(2*N)*r);
% 
% 
% lowerbound = sqrt(B(x, x, N));
% upperbound = 1000;
% temp = integral(f, lowerbound, upperbound);
% 
% 
% res = res-  temp;
% 
% end

function res = E(x, v, N)


V = B(v, v, N);
X = B(x, x, N);

res = sqrt(V);





temp = (pi/2 - atan(sqrt(X))) / sqrt(2*N);


res = res-  temp;

end
