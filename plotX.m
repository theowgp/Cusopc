function  plotX( sol, t, N, d )
%% PLOT X
% figure
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    YX(k) =  B(x, x, N);
end
plot(t, YX);
hold all
title('X(t)');
end

