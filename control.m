function res = control(x, v, i, N, d, R, A)


temp = zeros(1, d);
for j = 1:N
    if A(i, j)
%         temp = temp-    dV(norm(x(i, :) - x(j, :)), R) * (x(i, :) - x(j, :)) / norm(x(i, :) - x(j, :));
        
        rij = norm(x(i, :) - x(j, :));
        if rij >  R/2
            temp = temp+    -2 * V(rij, R) * (v(i, :) - v(j, :));
        else
            temp = temp+     2 * V(rij, R) * (v(i, :) - v(j, :));
        end
    end
end

res = temp;

end

