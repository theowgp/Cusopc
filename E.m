function res = E(x, v, N)

V = B(v, v, N);
X = B(x, x, N);

res = sqrt(V);


temp = (1/sqrt(2*N)) * (pi/2 - atan(sqrt(2*N*X)));

res = res-  temp;


end

