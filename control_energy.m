function res = control_energy(solx, dynamics, mesh, N, d, key)

temp = 0;
for k = 1:mesh.n+1
    [x, v] = convert(solx(:, k), N, d);
    for i = 1:N
        temp = temp+    norm(dynamics.control(x, v, i, key, 1))^2;
    end
end

res = temp * mesh.h;
end