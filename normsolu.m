function res = normsolu(solu, mesh)

res = sqrt(spsolu(solu, solu, mesh));

end
