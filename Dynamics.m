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
        

        
        
        
        function res = fx(obj, v)
            res = v;
        end
        
            
        
        function [res, alpha] = fv(obj, x, v, key)
            res = zeros(obj.N, obj.d);
            alpha = zeros(obj.N, 1);

            for i=1:obj.N
                temp = zeros(1, obj.d);
                for j=1:obj.N
                    temp = temp+  obj.a(norm(x(i, :) - x(j, :))) * (v(j, :) - v(i, :));
                end 
                if key == 1
%                     my control
                    alpha(i) = obj.DetermineAlpha(x, v, i);
                    res(i, :) = temp/obj.N + obj.control(x, v, alpha(i), i);
                else
                    if key == -1
%                         position averaging control
                        res(i, :) = temp/obj.N + obj.M*(obj.amean(x, x, i) - x(i, :));
                    else
%                         BFK control
                        res(i, :) = temp/obj.N + obj.M*(obj.amean(x, v, i) - v(i, :));
                    end
                    
                end
            end
        end
        
        function res = control(obj, x, v, alphai, i)
            res = alphai*(obj.amean(x, x, i) - x(i, :)) + (1 - alphai)*(obj.amean(x, v, i) - v(i, :));
            res = res * obj.M;
        end
        
        
        
        function res = DetermineAlpha(obj, x, v, i)
            t1 = norm(v(i, :) - obj.mean(v));
            t2 = norm(obj.amean(x, v, i) - obj.mean(v));
%             t3 = norm(x(i, :) - obj.mean(x));
            t4 = norm(obj.amean(x, x, i) - obj.mean(x));
            t5 = (x(i, :) - obj.mean(x))*(v(i, :) - obj.mean(v))';
            
            
            if t4*t1<t5 && t2>=t1
                res = 1;
            else
                if t4*t1>=t5 && t2<t1
                    res = 0;
                else
                    if t4*t1<t5 && t2<t1
                        if t4*t1 - t5 < t2 - t1 
                            res = 1;
                        else
                            res = 0;
                        end
                    else
                        if t4*t1>=t5 && t2>=t1
                           if t4*t1 - t5 < t2 - t1 
                                res = 1;
                           else
                                res = 0;
                           end 
                        end
                    end
                end
            end
                            
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
        

               
        
        
      
        
        
   
        
        
        function [res, alpha] = F(obj, argx, u)
            [x, v] = convert(argx, obj.N, obj.d);
            [fv, alpha] = obj.fv(x, v, u);
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(fv', [obj.N*obj.d, 1])];
        end
        
        
      
        
       
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

