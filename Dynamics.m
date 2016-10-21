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
        
            
        
        function res = fv(obj, x, v, key, ckey)
            res = zeros(obj.N, obj.d);
            

            for i=1:obj.N
                temp = zeros(1, obj.d);
                for j=1:obj.N
                    temp = temp+  obj.a(norm(x(i, :) - x(j, :))) * (v(j, :) - v(i, :));
                end 
                
%                 disp([i   v(i,:)]);
                res(i, :) = temp/obj.N + obj.control(x, v, i, key, ckey);
            end
        end
        
        function res = control(obj, x, v, i, key, ckey)
            if strcmp(key, 'my')
%                 global
%                 xhat = obj.mean(x) - x(i, :);
%                 local
                xhat = obj.amean(x, x, i) - x(i, :);
                if norm(xhat) ~= 0 && ckey == 1
                    vi = v(i, :);
                    if isnan(vi)
                        warning('liat');
                    end
                    teta = acos(vi * xhat'/(norm(vi)*norm(xhat)));
                    phi = teta - pi/2;

                    eibar = obj.Eibar(vi, xhat);
                    
                    
                    
                    res = eibar * cos(phi) *  obj.Sigma(norm(xhat)) * obj.Eta(teta);
%                     res = eibar * cos(phi) * obj.Sigma(norm(xhat));
                    
                    if isnan(res)
                        warning('liat');
                    end
                else
                    res = 0;
                end



% %                 global
% %               res = obj.mean(x) - x(i, :);
% %                 local    
                
%                 res = obj.amean(x, x, i) - x(i, :);
            end
            
            if strcmp(key, 'BFK')
%                 global 
%                 res = obj.mean(v) - v(i, :);

%                 local
                res = obj.amean(x, v, i) - v(i, :);
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
            res = 1/(1 + r^2)^obj.delta;
%             res = obj.cutoff(r);
        end
        

        
        
               
        
        function res = F(obj, argx, u, ckey)
            [x, v] = convert(argx, obj.N, obj.d);
            fv = obj.fv(x, v, u, ckey);
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(fv', [obj.N*obj.d, 1])];
        end
        
        
        
        
        
%         normilized orthogonal to velocity vector pointing to the similar direction of xhat
        function res = Eibar(obj, vi, xhat)
            ei = [vi(2)  -vi(1)];
            temp = ei*xhat';
            
            if temp > 0
                res = ei/norm(ei);
            else
                if temp < 0
                    res = -ei/norm(ei);
                else
                    res = zeros(size(ei));
                end
            end
        end
        
          
        
        function res = Sigma(obj, r) % Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            res = 10*r;% Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             res = 3*r;% Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             res = r;
        end
        
        
        function res = Eta(obj, teta) 
            res = (1 - cos(teta))*0.5;
        end
        
       
        
              
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

