function res = ComputeBetaK(g, gnext, drct, mesh)

y = gnext - g;

sigma = y - 2*drct*normsolu(y, mesh)^2/spsolu(y, drct, mesh);

res = spsolu(sigma, gnext, mesh)/spsolu(y, drct, mesh);

end