function tof_map = threshold_picker(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%varargin{1}: rfdata
%varargin{2}:kgrid
%varargin{3}:frequency [kHz]
%‚µ‚«‚¢’l–@‚Å”gŒ`æ“ª‚ð’Šo‚·‚éD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%“‰Šú‰»““““““““““““““““““““““““
if nargin == 2 % ‘—M”g’†SŽü”g”F‚QMH‚š
    rfdata = varargin{1};
    kgrid = varargin{2};
elseif nargin == 3 %‘—M”g’†SŽü”g”F2 MHzˆÈŠO
    rfdata = varargin{1};
    kgrid = varargin{2};
    frq = varargin{3};
else
    error('Number of variable is invalid.');
end
dt = kgrid.dt;
[~,r_num,t_num] = size(rfdata);
%“ŒvŽZ““““““““““““““““““““““““
if nargin == 2 % ‘—M”g’†SŽü”g”F‚QMH‚š
    tof_map = calc_threshold(rfdata, dt, t_num, r_num,43);%2018/08/12(realistic_model)
    %     tof_map = calc_threshold(rfdata, dt, t_num, r_num,41);%2018/08/08(pattern_model)
elseif nargin == 3 %‘—M”g’†SŽü”g”F2 MHzˆÈŠO
    switch frq
        case 50
            tof_map = calc_threshold(rfdata, dt, t_num, r_num,1702+2/3);%2018/08/12(realistic_model)
        case 100
            tof_map = calc_threshold(rfdata, dt, t_num, r_num,836+2/3);%2018/08/12(realistic_model)
        case 200
            tof_map = calc_threshold(rfdata, dt, t_num, r_num,415+2/3);%2018/08/12(realistic_model)
        case 500
            tof_map = calc_threshold(rfdata, dt, t_num, r_num,165);%2018/08/12(realistic_model)
        case 1000
            tof_map = calc_threshold(rfdata, dt, t_num, r_num,82+2/3);%2018/08/12(realistic_model)
        otherwise
            error('Designated Freqency is invalid.')
    end
end
end

function tof_map = calc_threshold(rfdata, dt, t_num, r_num, offset)
rfdata_abs = abs(rfdata);
threshold_map = max(rfdata_abs) * 10e-02;%2018/08/12(realistic_model)
threshold_map = reshape(threshold_map,r_num,t_num);
tof_map = zeros(r_num,t_num);
for transmit_i = 1:t_num
    for receive_j = 1:r_num
        ind_upper_threshold = (rfdata_abs(:,receive_j,transmit_i)>threshold_map(receive_j,transmit_i));
        k = find(ind_upper_threshold, 1 ) - offset;%2018/08/12(realistic_model)
        tof_map(receive_j,transmit_i) = k *dt;
    end
end
end