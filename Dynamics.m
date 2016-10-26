classdef Dynamics
       
    properties
       N;
       d;
        
       gamma;
       
       delta;
       
      
       
       eps;
       
       R;
       hysteresis;
       
       M;
       
       cp;
       
       A
        
    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, gamma, delta,  M, A, R, hysteresis)
            obj.N = N;
            obj.d = d;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.R = R;
            obj.hysteresis = hysteresis;
            obj.M = M;
            
            obj.A = A;
         
            
            %cutoff precision
%             obj.cp = 0.0000001;% worse than 0
            obj.cp = 0;
        end
        

        
        
        
        function res = fx(obj, v)
            res = v;
        end
        
            
        
        function [obj, res] = fv(obj, x, v, key)
            res = zeros(obj.N, obj.d);
            
            for i=1:obj.N
                [obj, ctrl] = obj.control(x, v, i, key);
                res(i, :) = ctrl;
            end
        end
        
        
        
        
        
        function [obj, res] = control(obj, x, v, i, key)
            if strcmp(key, 'my')
                    res = obj.control_my(x, v, i);
                else
                    if strcmp(key, 'BFK')
                        res = obj.control_BFK(x, v, i);
                    else
                        if strcmp(key, 'ZJP')
                            [obj, res] = obj.control_ZJP(x, v, i);
                        end
                    end
            end
        end
        
        
        
        function [obj, res] = control_ZJP(obj, x, v, i)
            res = zeros(1, obj.d);
            obj = obj.update_A(x);
            obj.A;
            
            for j = 1:obj.N
                if i ~= j && obj.A(i, j) == 1
                    res = res-  obj.A(i, j) * (v(i, :) - v(j, :));
                    res = res-  obj.A(i, j) * obj.GVijxi(x, i, j);
                end
            end
        end
        
        function res = GVijxi(obj, x, i, j)
%             disp([norm(x(1, :) - x(2, :))  norm(x(1, :) - x(3, :))  norm(x(2, :) - x(3, :))]);
            res =      -2 * x(i, :) / norm(x(i, :) - x(j, :))^4;
            
            res = res+  2 * (x(i, :) - x(j, :)) / (obj.R^2 - norm(x(i, :) - x(j, :))^2)^2;
        end
        
        function obj = update_A(obj, x)
            for i = 1:obj.N
                for j = 1:obj.N
                    if i ~= j
                        if norm(x(i, :) - x(j, :)) < obj.hysteresis
                            obj.A(i, j) = 1;
                        else
                            if norm(x(i, :) - x(j, :)) >= obj.R
                                obj.A(i, j) = 0;
                            end
                        end
                    end
                end
            end
        end
        
        
        
        
        function res = control_my(obj, x, v, i)
%                 global
%                 xhat = obj.mean(x) - x(i, :);
%                 local
                xhat = obj.amean(x, x, i) - x(i, :);
                if norm(xhat) ~= 0
                    vi = v(i, :);
                    
% %                     calculate angle between xhat and vi
%                     if isnan(vi)
%                         warning('liat');
%                     end
%                     teta = acos(vi * xhat'/(norm(vi)*norm(xhat)));
                    
%                     the direction of the cange of velocity
                    res = xhat - vi;

% %                     if it is nonzero then normilize it
%                     if norm(res) == 0
%                         warning('liat');
%                     else
%                         res = res/norm(res);
%                     end

% %                     multiply by force factors
%                     res = res * obj.Sigma(norm(xhat)) * obj.Eta(teta); 
                else
                    res = zeros(1, obj.d);
                end
          
                res = res * obj.M;
        end
        
        
        
        
        function res = control_BFK(obj, x, v, i)
%                 global 
%                 res = obj.mean(v) - v(i, :);

%                 local
                res = obj.amean(x, v, i) - v(i, :);
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
            res = 1/(1 + r^2)^obj.delta;
%             res = obj.cutoff(r);
        end
        
        
        
                

        
        
               
        
        function [obj, res] = F(obj, argx, key)
            [x, v] = convert(argx, obj.N, obj.d);
            [obj, fv] = obj.fv(x, v, key);
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(fv', [obj.N*obj.d, 1])];
        end
        
        
           
    end 
    
    
    methods(Static)
        
        
%         force factor depending on the distance between agents
        function res = Sigma(r) 
%             res = 10*r;% Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             res = 3*r;% Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            res = r;
%               res = obj.a(r);% bad idea!!!!
        end
 
        
%         force factor depending on how different the velocity vector is from the direction in which it shoud be changed
        function res = Eta(teta) 
            res = (1 - cos(teta))*0.5;
        end
        
        
    end
    
end

