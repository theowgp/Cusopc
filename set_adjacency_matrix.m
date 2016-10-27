function A = set_adjacency_matrix(x, N, R)
            A = zeros(N);
            for i = 1:N
                for j = 1:N
                    if i ~= j
                        if norm(x(i, :) - x(j, :)) < R
                            A(i, j) = 1;
                        end
                    end
                end
            end

end

