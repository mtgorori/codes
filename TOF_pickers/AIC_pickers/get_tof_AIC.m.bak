%%%% Load kgrid and sensor structure data in advance. %%%%
function tof_map = get_tof_AIC(varargin)
% varargin{1} is the receive data of objective
% varargin{2} is the receive data of water
% varargin{3} is the kgrid data
% varargin{4} is the sensor data
% varargin{5} is the half-size of window
display('start running AIC method for TOF');
tic;
tof_map = tof_AIC(varargin{1}, varargin{3}, varargin{4},varargin{5}) - tof_AIC(varargin{2}, varargin{3}, varargin{4},varargin{5});
display(['AIC method has gotten TOF diff and used ',num2str(toc),' s']);
end



function tof_map = tof_AIC(rfdata,kgrid,sensor,window_size_half)
v_water = 1540;
dt = kgrid.dt;
[~,t_num] = size(sensor.mask);
[frame_num,~,~] =size(rfdata);
element_dis = zeros(t_num,t_num);

for ii = 1:t_num % transmit
    for jj = 1: t_num % receiver
        if ii == jj
            continue
        end
        element_dis(jj,ii) = sqrt((sensor.mask(1,ii)-sensor.mask(1,jj))^2 + (sensor.mask(2,ii)-sensor.mask(2,jj))^2);% element distance to the first elements
    end
end

tof_water = element_dis ./ v_water;
% window_size_half = 75;
tof_map = zeros(t_num,t_num);

for transmit_i = 1:t_num
    curr_rfdata = rfdata(:,:,transmit_i);
    curr_rfdata_extent = cat(1,zeros(size(curr_rfdata((frame_num -window_size_half+1):end,:))),curr_rfdata);
    window_center = round(tof_water(:,transmit_i)/dt) + window_size_half;
    window_start = window_center - window_size_half +1;
    window_end = window_center + window_size_half;
    
    for receive_j = 1:t_num
        if abs(transmit_i - receive_j) <=3 
            continue
        end
        scan_region = curr_rfdata_extent(window_start(receive_j):window_end(receive_j),receive_j);
        ind = aic_pick(scan_region);
        tof_map(receive_j,transmit_i) = ((ind-1) + window_start(receive_j) - window_size_half) * dt;
    end
end
tof_map(tof_map<0) = 0;
end