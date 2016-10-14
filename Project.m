% function res = Project(solu, step, drct, N)
% [Nu, n, s] = size(solu);
% 
% res = solu +step*drct;
% 
% 
% 
% for i=1:N
%     for k=1:n
%         for l=1:s
%             if sqrt(res(2*i-1, k, l)^2 + res(2*i, k, l)^2) > 1
%                 res(2*i-1, k, l) = res(2*i-1, k, l)/sqrt(res(2*i-1, k, l)^2 + res(2*i, k, l)^2);
%                 res(2*i, k, l) = res(2*i, k, l)/sqrt(res(2*i-1, k, l)^2 + res(2*i, k, l)^2);
%             end
%             
%             if res(2*i, k, l)<0
%                 res(2*i, k, l) = 0;
%             end
%             
%             if res(2*i-1, k, l)<0
%                 res(2*i-1, k, l) = 0;
%             end
%             
%         end
%             
%     end
% end
% 
% 
% end
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