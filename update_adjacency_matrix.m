function A = update_adjacency_matrix(A, x, N, R, Rh)
            
            for i = 1:N
                for j = 1:N
                    if i ~= j
                        if norm(x(i, :) - x(j, :)) < Rh
                            A(i, j) = 1;
                        else
                            if norm(x(i, :) - x(j, :)) >= R
                                A(i, j) = 0;
                            end
                        end
                    end
                end
            end

end

