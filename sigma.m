function res = sigma(x, v, i, R, N, d)

if (x(i, :) - meani(x, i, R, N, d)) * (v(i, :) - meani(v, i, R, N, d))' > 0
    res = 1;
else
    res = 0;
end



end

