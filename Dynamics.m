classdef Dynamics
       
    properties
       N;
       d;
        
       gamma;
       
       delta;
       
       eps;
       
       R;
       
       M;
        
    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, gamma, delta, M, R)
            obj.N = N;
            obj.d = d;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.R = R;
            obj.M = M;
        end
        

        
        
        
        function res = fx(obj, v)
            res = v;
        end
        
            
        
        function res = fv(obj, x, v, u)
%         function res = fv(obj, x, v)
            res = zeros(obj.N, obj.d);

            for i=1:obj.N
                temp = zeros(1, obj.d);
                for j=1:obj.N
                    temp = temp+  obj.a(norm(x(i, :) - x(j, :))) * (v(j, :) - v(i, :));
                end 
                temp = temp/obj.N + obj.control(x, v, u, i);
                res(i, :) = temp;
            end
        end
        
        function res = control(obj, x, v, u, i)
            temp = u(i)*(obj.amean(x, x, i) - x(i, :)) + (1 - u(i))*(obj.amean(x, v, i) - v(i, :));
            nfactor = obj.M /norm(temp);
            res = temp * nfactor;
        end
        
        function res = amean(obj, x, w, i)
            temp = zeros(1, obj.d);
            for j=1:obj.N
                psiij = obj.cutoff(norm(x(i, :) - x(j, :)));
                temp = temp+    psiij*w(j, :);
            end
            
            res = temp;
            nfactor = obj.numofn(x, i);
            if nfactor > 0
                res = res/nfactor;
            end
        end
        
        function res = numofn(obj, x, i)
            res = 0;
            for j=1:obj.N
                psiij = obj.cutoff(norm(x(i, :) - x(j, :)));
 %                 if(psiij ~= 0)
                if(psiij > 0.00000001)
                    res = res+  1; 
                end
            end
        end
        
        % non normilized cutoff function
        function res = cutoff(obj, r)
            res = 0;
            if r < obj.R
                res = exp(-obj.R^2/(obj.R^2-r^2));
            end
        end
        
        % non normilized derivative of the cutoff function
        function res = dcutoff(obj, r)
            res = 0;
            if r < obj.R
                res = -2 * r * obj.R^2 * exp(-obj.R^2/(obj.R^2-r^2)) / (obj.R^2-r^2)^2;
            end
        end
        
        function res = a(obj, r)
            res = 1/(1 + r^2)^obj.delta;
        end
        
        function res = da(obj,  r)
            res = -obj.delta*2*r / (1 + r^2)^(1 + obj.delta);
        end
        
               
        
        
        function res = fz(obj, v)
            res = B(v, v, obj.N);
        end
        
        
        
        
        
        
        function res = dfvdx(obj, x, v, u, i, k)
            temp = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                temp = temp+    (v(j, :) - v(i, :))' * obj.da(norm(x(i, :) - x(j, :))) * obj.dnorm(x, i, j, k);
            end
            temp = temp/obj.N;
            
            res = temp + obj.dcontroldx(x, v, u, i, k);
        end
        
        function res = dnorm(obj, x, i, j, k)
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
        
        function res = dcontroldx(obj, x, v, u, i, k)
            res = u(i)*(obj.dameanxdx(x, i, k) - obj.dxdx(i, k)) + (1 - u(i))*obj.dameanvdx(x, v, i, k);
            temp = u(i)*(obj.amean(x, x, i) - x(i, :)) + (1 - u(i))*(obj.amean(x, v, i) - v(i, :));
            nfactor = obj.M /norm(temp);
            res = res * nfactor;
        end
        
        function res = dameanxdx(obj, x, i, k)
            temp = zeros(obj.d);
            for j=1:obj.N
                temp = temp+    x(j, :)'*obj.dcutoff(norm(x(i, :) - x(j, :))) * obj.dnorm(x, i, j, k)   +   obj.cutoff(norm(x(i, :) - x(j, :))) * obj.dxdx(j, k);
            end
            
            res = temp;
            nfactor = obj.numofn(x, i);
            if nfactor > 0
                res = res/nfactor;
            end
        end
        
        function res = dameanvdx(obj, x, v, i, k)
            temp = zeros(obj.d);
            for j=1:obj.N
                temp = temp+    v(j, :)'*obj.dcutoff(norm(x(i, :) - x(j, :))) * obj.dnorm(x, i, j, k);
            end
            
            res = temp;
            nfactor = obj.numofn(x, i);
            if nfactor > 0
                res = res/nfactor;
            end 
        end
        
        function res = dxdx(obj, i, k)
            res = zeros(obj.d);
            if(k == i)
                res = eye(obj.d);
            end
        end
        
        
        
        
        
        
        
        function res = dfvdv(obj, x, v, u, i, k)
            temp = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                temp = temp+   obj.a(norm(x(i, :) - x(j, :))) * (obj.dxdx(j, k) - obj.dxdx(i, k));
            end
            temp = temp/obj.N;
            
            res = temp + obj.dcontroldv(x, v, u, i, k);
        end
        
        function res = dcontroldv(obj, x, v, u, i, k)
            res = (1 - u(i)) * (obj.dameanvdv(x, v, i, k) - obj.dxdx(i, k));
            temp = u(i)*(obj.amean(x, x, i) - x(i, :)) + (1 - u(i))*(obj.amean(x, v, i) - v(i, :));
            nfactor = obj.M /norm(temp);
            res = res * nfactor;
        end
        
        
        function res = dameanvdv(obj, x, v, i, k)
            temp = zeros(obj.d);
            for j=1:obj.N
                temp = temp+    obj.cutoff(norm(x(i, :) - x(j, :))) * obj.dxdx(j, k);
            end
            
            res = temp;
            nfactor = obj.numofn(x, i);
            if nfactor > 0
                res = res/nfactor;
            end
        end
            
        
        
        
        
        
        
        
        
        
        
        
        function res = dfvdu(obj, x, v, u, i, k)
            res = zeros(obj.d, 1);
            if k == i
                temp1 = u(i)*(obj.amean(x, x, i) - x(i, :)) + (1 - u(i))*(obj.amean(x, v, i) - v(i, :));
                temp2 = obj.amean(x, x, i) - x(i, :) - obj.amean(x, v, i) + v(i, :);
                                
                res = temp2/norm(temp1)^2 - temp1*(temp1*temp2')/norm(temp1)^3;
                res = res';
            end
        end
        
        
        
        
        function [x, v, z, u] = convert(obj, argx, argu)
            x = reshape(argx(1 : obj.N*obj.d), [obj.d, obj.N])';
            v = reshape(argx(obj.N*obj.d+1 : 2*obj.N*obj.d), [obj.d, obj.N])';
            z = argx(2*obj.N*obj.d + 1);
            u = reshape(argu, [1, obj.N])';
        end
    
        
        
        function res = F(obj, argx, argu)
            [x, v, z, u] = obj.convert(argx, argu);
            
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(obj.fv(x, v, u)', [obj.N*obj.d, 1]);     obj.fz(v)];
        end
        
        
        function res = GuF(obj, argx, argu)
            [x, v, z, u] = obj.convert(argx, argu);
            
            N = obj.N;
            d = obj.d;
            
            res = zeros(2*N*d + 1, N);
            
            temp = zeros(N*d, N);
            
            %dfvdu
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, k) = obj.dfvdu(x, v, u, i, k); 
                end
            end
            res(N*d+1:2*N*d, 1:N) = temp;
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
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdx(x, v, u, i, k); 
                end
            end
            res(N*d+1:2*N*d, 1:N*d) = temp;
            
            %dfvdv
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdv(x, v, u, i, k);
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
        
       
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

