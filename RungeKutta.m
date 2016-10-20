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
        
        
        
        function [solx, soly] = solve_forward_equation(obj, key, ckey)
            solx = zeros(obj.N, obj.grid.n+1);
            solx(:, 1) = obj.arg0;
            soly = zeros(obj.N, obj.grid.n, obj.s);

            
            for k=1:obj.grid.n
                solx(:, k+1) = solx(:, k);
                for i=1:obj.s
                   soly(:, k, i) = solx(:, k);
                   for j=1:obj.s
                       if i>j
                           soly(:, k, i) = soly(:, k, i) + obj.grid.h * obj.A(i, j) * obj.dynamics.F(soly(:, k, j), key , ckey);
                       end
                   end
                   solx(:, k+1) = solx(:, k+1) + obj.grid.h*obj.b(i) * obj.dynamics.F(soly(:, k, i), key, ckey);
                   [x, v] = convert(solx(:, k+1), obj.N/2/2, 2);
                   temp = E(x, v, obj.N/2/2);
                   if temp < 0
                       ckey = 0;
                   end
                end
            end
        end
        
        
        
        
        
        
    
        
          
        
    end
    
end