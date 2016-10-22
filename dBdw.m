function res = dBdw(w, k, N, d)
            temp = zeros(1, d);
            
            for i = 1:N
                temp = temp+    norm(w(k, :) - w(i, :)) * (w(k, :) - w(i, :));
            end
            
            res = 2*temp/N^2;
end
