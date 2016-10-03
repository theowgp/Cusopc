classdef Solver
        
    properties
        mesh;
        
        A;
        b;
        c;
        
        s;
        
        N;
        
        
    end
    
    methods
        
        
        function obj = Solver(A, b, c, s, N, T, n)
                        
            obj.A = A;
            obj.b = b;
            obj.c = c;
            obj.s = s;
            obj.N = N;            
            
            obj.mesh = Mesh(T, n);   

        end
        
     
        
        function [t, sol] = rk(obj, x0, f)
            sol = zeros(obj.mesh.n+1, obj.N);
            sol(1, :) = x0;
            soly = zeros(obj.mesh.n, obj.s, obj.N);

            for k=1:obj.mesh.n
                sol(k+1, :) = sol(k, :);
                for i=1:obj.s
                   soly(k, i, :) = sol(k, :);
                   for j=1:obj.s
                       if i>j
                           soly(k, i, :) = reshape(soly(k, i, :), [1, obj.N]) + obj.mesh.h*obj.A(i, j)*f(soly(k, j, :));
                       end
                   end
                       sol(k+1, :) = sol(k+1, :) + obj.mesh.h*obj.b(i)*f(soly(k, i, :));
                end
            end
            
            t = obj.mesh.t;
        end
        
        
        
        
        
        
        function [t, sol] = euler(obj, x0, f)
            sol = zeros(obj.mesh.n+1, obj.N);
            sol(1, :) = x0;

            for k=1:obj.mesh.n
                sol(k+1, :) = sol(k, :) + obj.mesh.h*f(sol(k, :));
            end
            
            t = obj.mesh.t;
        end


        
        
        
        

          
    end
    
    
    
    
end

