classdef Objective
    
    
    properties
        N;
        d
        
        alpha1;
        alpha2;
        alpha3;
    end
    
    methods
        function obj = Objective(N, d, alpha1, alpha2, alpha3)
            obj.N = N;
            obj.d = d;
            
            obj.alpha1 = alpha1;
            obj.alpha2 = alpha2;
            obj.alpha3 = alpha3;
        end
        
        
%         function res = phi(obj, z)
%             res = z;
%         end
        function res = phi(obj, argx)
            v = reshape(argx(obj.N*obj.d+1 : 2*obj.N*obj.d), [obj.d, obj.N])';
            z = argx(2*obj.N*obj.d + 1);
            res = z;
            res = res+  obj.alpha2*B(v, v, obj.N);
        end
         
        function res = Gxphi(obj, argx)
            v = reshape(argx(obj.N*obj.d+1 : 2*obj.N*obj.d), [obj.d, obj.N])';
            
            res = zeros(2*obj.N*obj.d + 1, 1);
            res(2*obj.N*obj.d + 1) = 1;
            
            temp = zeros(obj.N*obj.d, 1);
            for k = 1:obj.N
                temp((k-1)*obj.d+1:k*obj.d, 1) = dBdw(v, k, obj.N, obj.d); 
            end
            temp = temp*     obj.alpha2;
            
            res(obj.N*obj.d+1:2*obj.N*obj.d, 1) = temp;
        end
    end
    
end
