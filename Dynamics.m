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
        

        
        
        
        function res = fx(obj, v)
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
        
               
        
        
        
        function res = fz(obj, v)
            res = B(v, v, obj.N);
        end
        
        
        
        
        
        
        function res = dfvdu(obj, x, v, k, i)
            res = zeros(obj.d);
            if k == i
                temp = zeros(1, obj.d);
                for j = 1:obj.N
                    temp = temp+    obj.cutoff(norm(x(k, :) - x(j, :))) * (v(k, :) - v(j, :) + x(j, :) - x(k, :));
                end
                res = diag(temp);
            end
        end
        
        
        
        
        function [x, v, z, u] = convert(obj, argx, argu)
            x = reshape(argx(1 : obj.N*obj.d), [obj.d, obj.N])';
            v = reshape(argx(obj.N*obj.d+1 : 2*obj.N*obj.d), [obj.d, obj.N])';
            z = argx(2*obj.N*obj.d + 1);
            u = reshape(argu, [obj.d, obj.N])';
        end
    
        
        
        function res = F(obj, argx, argu)
            [x, v, z, u] = obj.convert(argx, argu);
            
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(obj.fv(x, v, u)', [obj.N*obj.d, 1]);     z];
        end
        
        
        function res = GuF(obj, argx, argu)
            [x, v, z, u] = obj.convert(argx, argu);
            
            N = obj.N;
            d = obj.d;
            
            res = zeros(2*N*d + 1, N*d);
            
            temp = zeros(N*d);
            
            %dfvdu
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdu(x, v, k, i); 
                end
            end
            res(N*d+1:2*N*d, 1:N*d) = temp;
        end
        
        
        
        function res = GxF(obj, argx, argu)
            [x, v, z, u] = obj.convert(argx, argu);
            
            N = obj.N;
            d = obj.d;
            
            temp = zeros(N*d);
            
            res = zeros(2*N*d + 1);
            
            %dfxdv
            res(1:N*d,    N*d+1:2*N*d) = eye(N*d);
                                   
            %dfvdx
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdx(x, v, u, k, i); 
                end
            end
            res(N*d+1:2*N*d, 1:N*d) = temp;
            
            %dfvdv
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdv(x, u, k, i);
                end
            end
            res(N*d+1:2*N*d, N*d+1:2*N*d) = temp;
            
            %dfzdv
            temp = zeros(1, N*d);
            for k = 1:N
                temp((k-1)*d+1:k*d) = obj.dfzdv(v, k); 
            end
            res(2*N*d+1, N*d+1:2*N*d) = temp;
        end
        
        
        
        
        
        
        
        
        
        
        
        
        function res = dfzdv(obj, v, k)
            temp = zeros(1, obj.d);
            
            for i = 1:obj.N
                temp = temp+    norm(v(k, :) - v(i, :)) * (v(k, :) - v(i, :));
            end
            
            res = 2*temp/obj.N^2;
        end
        
        
        
        
        
        
        function res = dfvdv(obj, x, u, k, i)
            tempm1 = zeros(obj.d, obj.d);
            tempm2 = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                aij = obj.a(x, i, j);
                tempm1 = tempm1+    aij * obj.dtemp(k, i, j);
                tempm2 = tempm2+    aij * (- obj.dtemp(k, i, j));                           
            end
            % in this loop I multiply control component corresponding to the term tempm2 component 
            for s = 1:obj.d
                tempm2(s, :) = tempm2(s, :) * u(i, s); 
            end
            res = 2*tempm1 + tempm2;
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
            if k == j
                res =  eye(obj.d);
            else
                if k == i
                    res = - eye(obj.d);
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
            if i ~= j
                if k == i
                    res = x(k, :)/norm(x(k, :) - x(j, :));
                else
                    if k == j
                        res = -x(k, :)/norm(x(i, :) - x(k, :));
                    end
                end
            end
        end
        
        function res = dcutoffdx(obj, x, k, i, j)
            res = obj.dcutoff(norm(x(i, :) - x(j, :))) * obj.dnorm(x, k, i, j);
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

