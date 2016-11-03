function res = control(x, i, N, d, R, A)


temp = zeros(1, d);
for j = 1:N
    if A(i, j)
        temp = temp-    dV(norm(x(i, :) - x(j, :)), R) * (x(i, :) - x(j, :)) / norm(x(i, :) - x(j, :));
    end
end

res = temp;

end

