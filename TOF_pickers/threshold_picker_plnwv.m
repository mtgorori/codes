function tof_map = threshold_picker_plnwv( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%varargin{1}: rfdata
%varargin{2}:kgrid
% varargin{3}:number_of_channels;'all','single'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%しきい値法で波形先頭を抽出する．
%初期化
rfdata = varargin{1};
kgrid = varargin{2};
dt = kgrid.dt;


if nargin == 2
    [~,r_num] = size(rfdata);
    tof_map = zeros(1,r_num);
    rfdata_abs = abs(rfdata);
    % threshold_map = max(rfdata_abs) * 200e-04;%2018/08/08(pattern_model)
    threshold_map = max(rfdata_abs) * 10e-02;%2018/08/12(realistic_model)
    % threshold_map = reshape(threshold_map,r_num,t_num);
    for receive_j = 1:r_num
        ind_upper_threshold = (rfdata_abs(:,receive_j)>threshold_map(1,receive_j));
        %         k = find(ind_upper_threshold, 1 ) - 41;%2018/08/08(pattern_model)
        %     k = find(ind_upper_threshold, 1 ) - 43;%2018/08/12(realistic_model)
        k = find(ind_upper_threshold, 1 ) - 211;%2018/08/23(realistic_model_plwv)
        tof_map(1,receive_j) = k *dt;
    end
elseif nargin == 3
    [~,~,r_num] = size(rfdata);
    switch varargin{3}
        case 'all'
            tof_map = zeros(1,r_num);
            rfdata_abs = abs(rfdata);
            % threshold_map = max(rfdata_abs) * 200e-04;%2018/08/08(pattern_model)
            threshold_map = max(rfdata_abs) * 10e-02;%2018/08/12(realistic_model)
            % threshold_map = reshape(threshold_map,r_num,t_num);
            for receive_j = 1:r_num
                ind_upper_threshold = (rfdata_abs(:,receive_j)>threshold_map(1,receive_j));
                %         k = find(ind_upper_threshold, 1 ) - 41;%2018/08/08(pattern_model)
                %     k = find(ind_upper_threshold, 1 ) - 43;%2018/08/12(realistic_model)
                k = find(ind_upper_threshold, 1 ) - 211;%2018/08/23(realistic_model_plwv)
                tof_map(1,receive_j) = k *dt;
            end
        case 'single'
            % rfdata (1chずつ送信)を合計して単一ch送受信を模擬する．
            tof_map = zeros(1,1);
            rfdata_sum = sum(rfdata,3);
            rfdata_sum = sum(rfdata_sum(:,101:200),2);
            rfdata_sum = rfdata_sum/(r_num^2);
            rfdata_abs = abs(rfdata_sum);
            thereshold = max(rfdata_abs) * 10e-02;
            ind_upper_threshold = rfdata_abs>thereshold;
            k = find(ind_upper_threshold,1) -170.7;
            tof_map(1,1) = k *dt;
    end
end
end