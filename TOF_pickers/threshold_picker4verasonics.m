function tof_map = threshold_picker4verasonics(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%varargin{1}: rfdata
%varargin{2}:dt
%varargin{3}:frequency [kHz]
%しきい値法で波形先頭を抽出する．
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%％初期化％％％％％％％％％％％％％％％％％％％％％％％％
rfdata = varargin{1};
dt = varargin{2};
frq = varargin{3};
[~,r_num,t_num] = size(rfdata);

%％計算％％％％％％％％％％％％％％％％％％％％％％％％
    switch frq
        case 5000
            tof_map = calc_threshold(rfdata, dt, t_num, r_num,3);
        otherwise
            error('Designated Freqency is invalid.')
    end
end

function tof_map = calc_threshold(rfdata, dt, t_num, r_num, offset)
rfdata_abs = abs(rfdata);
threshold_map = max(rfdata_abs)/50;
threshold_map = reshape(threshold_map,r_num,t_num);
tof_map = zeros(r_num,t_num/2);
for transmit_i = t_num/2+1:t_num
    for receive_j = 1:r_num
        ind_upper_threshold = (rfdata_abs(:,receive_j,transmit_i)>threshold_map(receive_j,transmit_i));
        k = find(ind_upper_threshold, 1 ) - offset;
        tof_map(receive_j,transmit_i) = k *dt;
    end
end
end