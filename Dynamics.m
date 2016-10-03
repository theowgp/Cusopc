classdef Dynamics
       
    properties
       N;
       d;
        
       gamma;
       
       delta;
       
       eps;
       
       R;
       
    
    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, gamma, delta, R)
            obj.N = N;
            obj.d = d;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.R = R;
        end
        

        
        
        
        function res = fx(obj, x, v)
            res = v;
        end
        
        
        
        
        function res = fv(obj, x, v)
            res = zeros(obj.N, obj.d);

            for i=1:obj.N
                for j=1:obj.N
                    res(i, :) = res(i, :)+  obj.a(norm(x(i, :) - x(j, :))) * (v(j, :) - v(i, :))     +    obj.u(i, x, v);
                end
                res(i, :) = res(i, :) / obj.N;
            end
            
        end
        
        
        
        function res = f(obj, x, v)
            res = [obj.fx(x, v); obj.fv(x, v)];
        end

        
        
    
        function res = a(obj, r)
            res = 1 / (1+r^2)^obj.delta;
        end
        
        
        function res = u(obj, i, x, v)
            res = -obj.gamma * (v(i, :) - obj.mean(i, x, v, obj.R)) ;
            res = res -obj.gamma * (x(i, :) - obj.mean(i, x, x, obj.R)) ;
        end
        
                
        function res = eta(obj, x, R)
            res = 0;
            max1 = 0;
            for i=1:obj.N
                max0 = 0;
                for k=1:obj.N
                    max0 = max0+ obj.khi(norm(x(i, :) - x(k, :)), R);
                end
                if max1<max0
                    max1 = max0;
                end
            end
            res = max1;
        end
        
        
        
        function res = mean(obj, i, x, v, R)
            res = zeros(1, obj.d);
            for j=1:obj.N
                res = res+ obj.khi(norm(x(i, :) - x(j, :)), R) * v(j, :);
            end
            res = res/obj.eta(x, R);
        end
    
        
    
    end 
    
    
    methods(Static)
        
        
        
        
        function res = khi(r, R)
            res = 0;
            if r<=R
                res = 1;
            end
        end
        
        
        
    end
    
end

