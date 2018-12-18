%%%% Load kgrid and sensor structure data in advance. %%%%
function tof_map = get_tof_AIC_from_singleRFData(varargin)
%背景RFデータを用いずにTTM(Travel Time Map, not TTDM(Tiravel Time Difference
%Map))を作成するための関数．[2018-04-27]
% varargin{1} is the receive data of objective
% varargin{2} is the kgrid data
% varargin{3} is the half-size of window//recommend:75
display('start running AIC method for TOF'); %#ok<DISPLAYPROG>
tic;
tof_map = tof_AIC(varargin{1}, varargin{2}, varargin{3});
display(['AIC method has gotten TOF and used ',num2str(toc),' s']); %#ok<DISPLAYPROG>
end



function tof_map = tof_AIC(rfdata,kgrid,window_size_half)
%初期化
dt = kgrid.dt;
[frame_num,r_num,t_num] = size(rfdata);
tof_map = zeros(r_num,t_num);

for transmit_i = 1:t_num
    %計算窓を全パターンで設定できるように，先頭にwindowhalfsize分だけ0を代入した．[2018-04-27]
    curr_rfdata = rfdata(:,:,transmit_i);
    curr_rfdata_extent = cat(1,zeros(size(curr_rfdata((frame_num-window_size_half+1):end,:))),curr_rfdata);
    
    for receive_j = 1:r_num
%         if abs(transmit_i - receive_j) <=3 
%             continue
%         end
        %ピークが複数出てきたときの対策．[2018-04-27]
        [pk_wave,window_center] = findpeaks(abs(hilbert(curr_rfdata_extent(:,receive_j))),'MinPeakHeight',0.04,'NPeaks',1);
        window_start = window_center - window_size_half +1;
        window_end = window_center;
        %ピークに比べてある程度小さな絶対値をとる frameで計算窓を打ち切る．[2018-04-27]
        while (abs(hilbert(curr_rfdata_extent(window_end,receive_j))) > (1/30)*pk_wave)
            window_end = window_end + 1;
        end
        scan_region = curr_rfdata_extent(window_start:window_end,receive_j);
        ind = aic_pick(scan_region);
        tof_map(receive_j,transmit_i) = ((ind-1) + window_start - window_size_half) * dt;
    end
end
tof_map(tof_map<0) = 0;
end