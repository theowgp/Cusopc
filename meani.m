function res = meani(v, i, R, N, d)

temp = zeros(1, d);
nfactor = 0;
for j = 1:N
    if norm(v(i, :) - v(j, :)) <= R
        temp = temp+    v(j, :);
        nfactor = nfactor+  1;
    end
end

res = temp/nfactor;
end

