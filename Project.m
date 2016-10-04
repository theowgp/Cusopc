function res = Project(solu, step, drct)
[N, n, s] = size(solu);

res = solu +step*drct;



for i=1:N
    for k=1:n
        for l=1:s
            if res(i, k, l)>1
                res(i, k, l) = 1;
            else
                if res(i, k, l)<0
                    res(i, k, l) = 0;
                end
            end
        end
            
    end
end


end
