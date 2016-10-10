function res = dBdv(v, k, N, d)
            temp = zeros(1, d);
            
            for i = 1:N
                temp = temp+    norm(v(k, :) - v(i, :)) * (v(k, :) - v(i, :));
            end
            
            res = 2*temp/N^2;
end