classdef Dynamics
       
    properties
       N;
       d;
        
       gamma;
       
       delta;
       
      
       
       eps;
       
       R;
       
       M;
       
       cp;
        
    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, gamma, delta,  M, R)
            obj.N = N;
            obj.d = d;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.R = R;
            obj.M = M;
         
            
            %cutoff precision
%             obj.cp = 0.0000001;% worse than 0
            obj.cp = 0;
        end
        

        
        
        
%         function res = fx(obj, v)
%             res = v;
%         end
        function res = fx(obj, x, v, key)
            res = zeros(obj.N, obj.d);
            
            if strcmp(key, 'my')
                for i=1:obj.N
                    res(i, :) = v(i, :) + obj.amean(x, x, i) - x(i, :);
                end
            end
            if strcmp(key, 'BFK')
                res = v;
            end
        end
        
            
        
        function res = fv(obj, x, v, key)
            res = zeros(obj.N, obj.d);
            

            for i=1:obj.N
                temp = zeros(1, obj.d);
                for j=1:obj.N
                    temp = temp+  obj.a(norm(x(i, :) - x(j, :))) * (v(j, :) - v(i, :));
                end 
%                 res(i, :) = temp/obj.N + obj.control(x, v, i, key);
%                 res(i, :) = temp/obj.N;
                res(i, :) = temp/obj.N + obj.amean(x, v, i) - v(i, :);
            end
        end
        
        function res = control(obj, x, v, i, key)
            X = B(x, x, obj.N);
            V = B(v, v, obj.N);
            gamma1 = obj.a(sqrt(2*obj.N*X))/sqrt(X);
            gamma2 = sqrt(V); 
              
            if strcmp(key, 'my')
%                 Et = E(x, v, obj.N);
%                 
%                 if( Et >= 0)
%                     res = (obj.mean(v) - v(i, :) - gamma1*x(i, :))*gamma2;
%                 else
%                     res = zeros(1, obj.d);
%                 end

%                 global
                res = (obj.mean(v) - v(i, :) - gamma1*x(i, :))*gamma2;
% 
%                 local
%                 res = (obj.amean(x, v, i) - v(i, :) - gamma1*x(i, :))*gamma2;
% 
%                 another
%                 res = (- v(i, :) - gamma1*x(i, :))*gamma2;

            end
            
            if strcmp(key, 'BFK')
%                 Et = E(x, v, obj.N)
                
%                 global 
                res = obj.mean(v) - v(i, :);

%                 local
%                 res = obj.amean(x, v, i) - v(i, :);
            end

            res = res * obj.M;
        end
        
        
        
               
        
        
        function res = amean(obj, x, w, i)
            temp = zeros(1, obj.d);
            for j = 1:obj.N
                psiij = obj.cutoff(norm(x(i, :) - x(j, :)));
                temp = temp+    psiij*w(j, :);
            end
            
            res = temp;
            nfactor = obj.numofn(x, i);
            res = res/nfactor;
        end
        function res = mean(obj, w)
            temp = zeros(1, obj.d);
            for j = 1:obj.N
                temp = temp+    w(j, :);
            end
            
            res = temp/obj.N;
        end
        
        function res = numofn(obj, x, i)
            res = 0;
            for j = 1:obj.N
                psiij = obj.cutoff(norm(x(i, :) - x(j, :)));
    %             if(psiij ~= 0)
                if(psiij > obj.cp)
                    res = res+  1; 
                end
            end
        end
        
        % normilized cutoff function
        function res = cutoff(obj, r)
            res = 0;
%             if r < obj.R
            if r < obj.R - obj.cp
                res = exp(-obj.R^2/(obj.R^2-r^2))   * exp(1);
            end
        end
        
   
        
        function res = a(obj, r)
            res = 1./(1 + r.^2).^obj.delta;
%             res = obj.cutoff(r);
        end
        

        
        
        
      
        
        
   
        
        
        function res = F(obj, argx, key)
            [x, v] = convert(argx, obj.N, obj.d);
            fv = obj.fv(x, v, key);
            res = [reshape(obj.fx(x, v, key)', [obj.N*obj.d, 1]);    reshape(fv', [obj.N*obj.d, 1])];
        end
        
        
      
        
       
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

