function  plotdX( sol, t, N, d )
%% PLOT X
% figure
for k = 1:length(t)
    [x, v] = convert(sol(k, :), N, d);
    YdX(k) =  2*B(x, v, N);
end
plot(t, YdX);
hold all
title('dX(t)');
end

