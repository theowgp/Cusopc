function res = F(argx, N, d, R, Rh)

global A

[x, v] = convert(argx, N, d);

A = update_adjacency_matrix(A, x, N, R, Rh);

tempx = v;

tempv = zeros(N, d);
for i = 1:N
    tempv(i, :) = zeros(1, d);
    for j = 1:N
        if A(i, j)
%             tempv(i, :) = tempv(i, :)-     (v(i, :) - v(j, :));
        end
    end
%     tempv(i, :) = tempv(i, :)+  control(x, v, i, N, d, R, A) * sigma(x, v, i, R, N, d);
    tempv(i, :) = tempv(i, :)+  control(x, v, i, N, d, R, A);
end

res = [reshape(tempx', [N*d, 1]); reshape(tempv', [N*d, 1])];

end