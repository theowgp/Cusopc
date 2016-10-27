function res = F(argx, N, d, R, Rh)

global A

[x, v] = convert(argx, N, d);

A = update_adjacency_matrix(A, x, N, R, Rh);

tempx = v;

tempv = zeros(N, d);

for i = 1:N
    tempv(i, :) = zeros(1, d);
    for j = 1:N
        if i ~= j && A(i, j) ~= 0
            tempv(i, :) = tempv(i, :)-     (v(i, :) - v(j, :));
            tempv(i, :) = tempv(i, :)-     dV(norm(x(i, :) - x(j, :)), R) * x(i, :) / norm(x(i, :) - x(j, :));
        end
    end
end

res = [reshape(tempx', [N*d, 1]); reshape(tempv', [N*d, 1])];

end