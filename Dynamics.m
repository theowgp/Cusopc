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
                        warning('sukabliat');
                    end
                    teta = acos(vi * xhat'/(norm(vi)*norm(xhat)));
                    ei = [v(i, 2), -v(i, 1)];


                    res = obj.Eibar(ei, xhat, teta, vi) * obj.Sigma(norm(x(i, :) - xhat)) * obj.Eta(teta);
                    if isnan(res)
                        warning('sukabliat');
                    end
                else
                    res = 0;
                end



% %                 global
% %                 res = obj.mean(x) - x(i, :);
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
            res = 1./(1 + r.^2).^obj.delta;
%             res = obj.cutoff(r);
        end
        

        
        
               
        
        function res = F(obj, argx, u, ckey)
            [x, v] = convert(argx, obj.N, obj.d);
            fv = obj.fv(x, v, u, ckey);
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(fv', [obj.N*obj.d, 1])];
        end
        
        
        
        
        
        
        
        function res = Eibar(obj, ei, xhat, teta, vi)
            temp = ei * xhat';
            
%             teta = acos(vi * xhat'/(norm(vi)*norm(xhat)));
            phi = obj.Phi(teta);
            
            if pi/2 < teta && teta <= pi
                if temp > 0
                    res = obj.Rotate(ei, phi, vi);
                else
                    if temp < 0
                        res = obj.Rotate(-ei, phi, vi);
                    else
                        res = -vi;
                    end
                end
            else
                if temp > 0
                    res = ei;
                else
                    if temp < 0
                        res = -ei;
                    else
                        res = -vi;
                    end
                end
            end
                
            
            if norm(res) == 0
                warning('sukabliat');
            end
            res = res/norm(res);
            if isnan(res)
                warning('sukabliat');
            end
        end
      
        
        
        
        
        function res = Sigma(obj, r) % Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            res = 3*r;% Achtung!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             res = r;
        end
        
        
        function res = Eta(obj, teta) 
            res = (1 - cos(teta))*0.5;
        end
        
        function res = CC(obj, phi)
            res = [cos(phi)  -sin(phi); sin(phi)  cos(phi)];
        end
        
        function res = Phi(obj, teta)
            res = 0;
            if pi/2 < teta && teta <= pi
                res = teta - pi/2;
            end
        end
        
        function res = Rotate(obj, ei, phi, vi)
            temp1 = obj.CC(phi)*ei';
            temp2 = obj.CC(-phi)*ei';
            
            if temp1' * vi' < 0
                res = temp1';
            else 
                res = temp2';
            end
        end
       
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

