classdef Mesh
        
    properties
         h;
         T;
         n;
         t;
         
    end
    
    methods
        function obj = Mesh(T, n)
            obj.T = T;
            obj.n = n;
            obj.h = T/n;
            obj.t = 0:obj.h:T;
        end
    
    end
    
end

