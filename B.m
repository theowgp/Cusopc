function res = B(x, y, N)

res = 0;
for i=1:N
    for j=1:N
        res = res+  (x(i, :) - x(j, :))*(y(i, :) - y(j, :))';
    end
end

res = res/(2*N^2);
end
