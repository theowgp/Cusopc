function res = dV(r, R)

res = -R^2/4/((r - R/2)^2  - R^2/4)^2 * 2 * (r - R/2);

% res = -2/r^3 + 2*r/(R^2 - r^2)^2;

% res = 2*(r - R/2);

end

