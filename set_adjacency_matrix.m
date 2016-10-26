function A = set_adjacency_matrix(x, N, R, hysteresis)
            A = zeros(N);
            for i = 1:N
                for j = 1:N
                    if i ~= j
                        if norm(x(i, :) - x(j, :)) < hysteresis
                            A(i, j) = 1;
                        end
                    end
                end
            end

end

