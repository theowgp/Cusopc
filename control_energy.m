function res = control_energy(t, sol, N, d, R, Rh)

A = zeros(N);

YE = zeros(length(t));

temp1 = 0;
for k =1:length(t)
    [x, v] = convert(sol(k, :), N, d);
    A = update_adjacency_matrix(A, x, N, R, Rh);
    
    temp2 = 0;
    for i = 1:N
        temp2 = temp2 + norm(control(x, v, i, N, d, R, A))^2;
    end
    temp1 = temp1 + temp2;
    
    YE(k) = temp2;
end

h = (t(end) - t(1)) / length(t);

res = temp1 * h;

% figure
% plot(t, YE);
% title('control energy');
end

