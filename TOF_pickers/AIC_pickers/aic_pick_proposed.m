function [ind] = aic_pick_proposed(x)

x = x - median(x); % remove median of x and window
% switch o
%     case {'to_peak'}
%         ind_peak = find(abs(x) == max(abs(x)));
%         xnew = x(1:ind_peak);
%     otherwise
%         xnew = x;
% end

junk = proposed_aicval(x);

if junk ~= 0
    ind = find(junk == min(junk)) + 1; % pick is one more than divide point
else
    ind = 0;
end
end


function [a] = proposed_aicval(x)
if ~isempty(x)
    n = length(x);
    a = zeros(n-1,1);
    for i=1:n-1
        %compute variance in first part
        s1 = 1/2 + log(std(x(1:i))) + (1/2)*log(2*pi);
        %compute variance in second part
        s2 =1 + log(std(x(i+1:n))) + (1/2)*log(2);
        a(i) = i*(s1) + (n-i+1)*(s2);
    end
else
    a = 0;
end
a(isinf(a)) = 0;
end