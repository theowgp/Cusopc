classdef Objective
    
    
    properties
        N;
        d
        
        dynamics;
                
    end
    
    methods
        function obj = Objective(N, d, dynamics)
            obj.N = N;
            obj.d = d;
            obj.dynamics = dynamics;
        end
        
        
%         function res = phi(obj, z)
%             res = z;
%         end
        function res = phi(obj, argx)
            v = reshape(argx(obj.N*obj.d+1 : 2*obj.N*obj.d), [obj.d, obj.N])';
            x = reshape(argx(1 : obj.N*obj.d), [obj.d, obj.N])';
            z = argx(2*obj.N*obj.d + 1);
            
            res = z;
            
            upperbound = 100;
            upperbound = upperbound/sqrt(2*obj.N);
            res = res+  E(x, v, obj.N, obj.dynamics, upperbound);  
        end
         
        function res = Gxphi(obj, argx)
            v = reshape(argx(obj.N*obj.d+1 : 2*obj.N*obj.d), [obj.d, obj.N])';
            x = reshape(argx(1 : obj.N*obj.d), [obj.d, obj.N])';
            
            XT = B(x, x, obj.N);
            VT = B(v, v, obj.N);
            
            res = zeros(2*obj.N*obj.d + 1, 1);
            
%             dphidx
            temp = zeros(obj.N*obj.d, 1);
            for k = 1:obj.N
                temp((k-1)*obj.d+1:k*obj.d, 1) = dBdv(x, k, obj.N, obj.d);  
            end
            res(obj.N*obj.d+1:2*obj.N*obj.d, 1) = 0.5 * temp * obj.dynamics.a(sqrt(2*obj.N*XT)) / sqrt(XT);
            res(1:obj.N*obj.d, 1) = temp;
            
%             dphidv
            temp = zeros(obj.N*obj.d, 1);
            for k = 1:obj.N
                temp((k-1)*obj.d+1:k*obj.d, 1) = dBdv(v, k, obj.N, obj.d);  
            end
            res(obj.N*obj.d+1:2*obj.N*obj.d, 1) = 0.5 * temp / sqrt(VT);
            
%             dphidz
            res(2*obj.N*obj.d + 1) = 1;
        end
    end
    
end
