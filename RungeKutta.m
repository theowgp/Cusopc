classdef RungeKutta
        
    properties
        grid;
        
        N;
        Nu;
        
        A;
        b;
        s;
        
        dynamics;  
       
                
        arg0;
     
    end
    
    methods
                
        function obj = RungeKutta(A, b, s, dynamics, arg0, N, Nu, T, n)
            obj.grid = Mesh(T, n);  
            
            obj.A = A;
            obj.b = b;
            obj.s = s;
            
            obj.N = N;
            obj.Nu = Nu;
            
            obj.dynamics = dynamics;
      
            obj.arg0 = arg0;
        end
        
        
        
        function [solx, soly] = solve_forward_equation(obj, key)
            solx = zeros(obj.N, obj.grid.n+1);
            solx(:, 1) = obj.arg0;
            soly = zeros(obj.N, obj.grid.n, obj.s);
            
            
            for k=1:obj.grid.n
%                 [x, v] = convert(solx(:, k), obj.N/4, 2);
%                 norm(v)
%                 Et = E(x, v, obj.N/4)
                
                solx(:, k+1) = solx(:, k);
                for i=1:obj.s
                   soly(:, k, i) = solx(:, k);
                   for j=1:obj.s
                       if i>j
                           f = obj.dynamics.F(soly(:, k, j), key);
                           soly(:, k, i) = soly(:, k, i) + obj.grid.h * obj.A(i, j) * f;
                       end
                   end
                   f = obj.dynamics.F(soly(:, k, j), key);
                   solx(:, k+1) = solx(:, k+1) + obj.grid.h*obj.b(i) * f;
                end
            end
        end
        
        
        
        
        
        
    
        
          
        
    end
    
end