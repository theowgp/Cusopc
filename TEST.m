


for k =1:length(t)
    [x, v] = convert(sol(k, :), N, d);
    A = update_adjacency_matrix(A, x, N, R, Rh);
    
    for i = 1:N
        xhat = meani(x, i, R, N, d) - x(i, :);
        test_value = control(x, i, N, d, R, A) * xhat'
    end
end