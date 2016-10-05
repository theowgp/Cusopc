function res = tempf(r, R)

for i = 1:length(r)
            res(i) = 0;
            if r(i) < R
                res(i) = exp(-R^2/(R^2-r(i)^2));
            end
end




end