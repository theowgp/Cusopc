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
        
            
        
%         function res = fv(obj, x, v, u)
        function res = fv(obj, x, v)
            res = zeros(obj.N, obj.d);

            for i=1:obj.N
                temp1 = zeros(1, obj.d);
                temp2 = zeros(1, obj.d);
                for j=1:obj.N
                    temp1 = temp1+  obj.a(x, i, j) * (v(j, :) - v(i, :));
                    temp2 = temp2+  obj.a(x, i, j) * (v(i, :) - v(j, :) + x(j, :) - x(i, :));
                end 
%                 res(i, :) = 2*temp1 + u(i)*temp2;
                res(i, :) = 2*temp1 + temp2;
            end
        end
        
        
    
        function res = f(obj, x, v)
            res = [obj.fx(x, v); obj.fv(x, v)];
        end
        
        
        
        
        % non normilized cutoff function 
        function res = cutoff(obj, r)
            res = 0;
            if r < obj.R
                res = exp(-obj.R^2/(obj.R^2-r^2));
            end
        end
        
        
        %coefficients aij in the system
        function res = a(obj, x, i,  j)
            res = obj.cutoff(norm(x(i, :) - x(j, :)));
            
            temp = 0;
            for s = 1:obj.N
                temp = temp+     obj.cutoff(norm(x(i, :) - x(s, :)));
            end
            
            res = res / temp;
        end
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

