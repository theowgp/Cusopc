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
        
            
        
        function res = fv(obj, x, v, u)
%         function res = fv(obj, x, v)
            res = zeros(obj.N, obj.d);

            for i=1:obj.N
                temp1 = zeros(1, obj.d);
                temp2 = zeros(1, obj.d);
                for j=1:obj.N
                    aij = obj.a(x, i, j);
                    temp1 = temp1+  aij * (v(j, :) - v(i, :));
                    temp2 = temp2+  aij * (v(i, :) - v(j, :) + x(j, :) - x(i, :));
                end 
                res(i, :) = 2*temp1 + u(i, :).*temp2;
%                 res(i, :) = 2*temp1 + temp2;
            end
        end
        
        
    
        function res = f(obj, x, v)
            res = [obj.fx(x, v); obj.fv(x, v)];
        end
        
        
        function res = fz(obj, x, v)
            res = B(v, v, obj.N);
        end
        
        
        
        function res = dfvdx(obj, x, v, u, k, i)
            tempm1 = zeros(obj.d, obj.d);
            tempm2 = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                tempm1 = tempm1+    (v(j, :) - v(i, :))' * obj.dadx(x, k, i, j);
                tempm2 = tempm2+    (v(i, :) - v(j, :) + x(j, :) - x(i, :))' * obj.dadx(x, k, i, j) + obj.dtemp(k, i, j) * obj.cutoff(norm(x(i, :) - x(j, :)));                           
            end
            % in this loop I multiply control component corresponding to the term tempm2 component 
            for s=1:obj.d
                tempm2(s, :) = tempm2(s, :) * u(i, s); 
            end
            res = 2*tempm1 + tempm2;    
        end
        
        function res = dtemp(obj, k, i, j)
            if k == i
                res = - eye(obj.d);
            else
                if k == j
                    res = eye(obj.d);
                else
                    res = zeros(obj.d);
                end
            end
        end
        
        function res = dadx(obj, x, k, i, j)
            temps = 0;
            tempv = zeros(1, obj.d);
            
            for s = 1:obj.N
                temps = temps+ obj.cutoff(norm(x(i, :) - x(s, :)));
                tempv = tempv+ obj.dcutoffdx(x, k, i, s);
            end
            
            res = obj.dcutoffdx(x, k, i, j) * temps - obj.cutoff(norm(x(i, :) - x(j, :))) * tempv;
            res = res / temps^2;
        end
        
        function res = dcutoff(obj, r)
            res = 0;
            if r < obj.R
                res = -2 * r * obj.R^2 * exp(-obj.R^2/(obj.R^2-r^2)) / (obj.R^2-r^2)^2;
            end
        end
        
        function res = dnorm(obj, x, k, i, j)
            res = zeros(1, obj.d);
            if k == i
                res = x(k, :)/norm(x(k, :) - x(j, :));
            else
                if k == j
                    res = -x(k, :)/norm(x(i, :) - x(k, :));
                end
            end
        end
        
        function res = dcutoffdx(obj, x, k, i, j)
            res = obj.dcutoff(norm(x(i, :) - x(j, :))) * obj.dnorm(obj, x, k, i, j);
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

