function [x, v] = convert(argx, N, d)
            x = reshape(argx(1 : N*d), [d, N])';
            v = reshape(argx(N*d+1 : 2*N*d), [d, N])';
            
           
end