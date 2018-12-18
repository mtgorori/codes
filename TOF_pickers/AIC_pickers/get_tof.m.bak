function [tof_map,tof_map2] = get_tof(x,kgrid,sensor)
%概略
%これは初期の関数．本来ならば計算窓の中心を与えてやらないといけない．ただ，これのいいところは，
%element distanceの出し方を参照できる点にある．それ以外は使い物にならない．[2018-04-27]

v_water = 1540;%水の音速[m/s]
% ring_radius = 50.e-3;% 半径
% rr =  ring_radius;
[~,t_num] = size(sensor.mask);%トランスデューサ数
% t = linspace(0, ((t_num-1)/t_num)*2*pi, t_num);%センサ位置角
% t_pos = zeros(2, t_num);%センサ位置
% t_pos(1,:) = rr*cos(t)+ring_cx;
% t_pos(2,:) = rr*sin(t)+ring_cy;
%     element_dis = zeros(t_num,t_num);
tof_map = zeros(t_num,t_num);
tof_map2 = zeros(t_num,t_num);
for ii = 1:t_num % transmit
    for jj = 1: t_num % receiver
        if ii == jj
            continue
        end
        element_dis = sqrt((sensor.mask(1,ii)-sensor.mask(1,jj))^2 + (sensor.mask(2,ii)-sensor.mask(2,jj))^2);% element distance to the first elements
        a = aic_pick(x(:,jj,ii),'to_peak') * kgrid.dt;
        tof_map(jj,ii) = a  - element_dis/v_water;
        tof_map2(jj,ii) = a ;
    end
end
end


function [ind] = aic_pick(x,o)

x = x - median(x); % remove median of x and window
switch o
    case {'to_peak'}
        ind_peak = find(abs(x) == max(abs(x)));
        xnew = x(1:ind_peak);
    otherwise
        xnew = x;
end

junk = aicval(xnew);

if junk ~= 0
    ind = find(junk == min(junk)) + 1; % pick is one more than divide point
else
    ind = 0;
end
end

function [a] = aicval(x)
if ~isempty(x)
    n = length(x);
    a = zeros(n-1,1);
    for i=1:n-1
        %compute variance in first part
        s1 = var(x(1:i));
        if s1 <= 0
            s1 = 0;
        else
            s1=log(s1);
        end
        %compute variance in second part
        s2 = var(x(i+1:n));
        if s2 <= 0
            s2 = 0;
        else
            s2=log(s2);
        end
        a(i) = i*(s1) + (n-i+1)*(s2);
    end
else
    a = 0;
end
end