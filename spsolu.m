function res = spsolu(solu, solv, mesh)
res = 0;

% for k = 1:mesh.n+1
for k = 1:mesh.n
    res = res+  solu(:, k, 1)'*solv(:, k, 1);
end

res = res*  mesh.h;
end

