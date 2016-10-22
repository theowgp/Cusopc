function [res, kA] = DetermineStepSize(rk, objective, mesh, solu, g, drct, sigma, limitA)

s = 10;

[solx, soly] = rk.solve_forward_equation(solu);
[solxA, solyA] = rk.solve_forward_equation(solu + s*drct);

kA = 0;
while objective.phi(solxA(:, end)) > objective.phi(solx(:, end))  + sigma * s * spsolu(drct, g, mesh)      &&     kA < limitA
    s = s * 0.5;
    [solxA, solyA] = rk.solve_forward_equation(solu + s*drct);
       
    kA = kA +1;
end




res = s;
end