classdef Objective
    
    
    properties
        N;
        d
    end
    
    methods
        function obj = Objective(N, d)
            obj.N = N;
            obj.d = d;
        end
        
        
%         function res = phi(obj, z)
%             res = z;
%         end
        function res = phi(obj, argx)
            z = argx(2*obj.N*obj.d + 1);
            res = z;
        end
        
        function res = Gxphi(obj)
            res = zeros(2*obj.N*obj.d + 1, 1);
            res(2*obj.N*obj.d + 1) = 1;
        end
    end
    
end
