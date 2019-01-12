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
tof_map = zeros(r_num/2,t_num/2);
    for transmit_i = 1:t_num/2
        for receive_j = 1:r_num/2
            [~,ind_max] = max(rfdata_abs(:,receive_j,transmit_i+t_num/2));
            ind_threshold = ind_max;
            while rfdata_abs(ind_threshold,receive_j,transmit_i+t_num/2) >= threshold_map(receive_j,transmit_i+t_num/2)
                ind_threshold = ind_threshold -1;
            end
            ind_threshold = ind_threshold - offset;
            tof_map(receive_j,transmit_i) = (ind_threshold - offset)* dt;
        end
    end
end