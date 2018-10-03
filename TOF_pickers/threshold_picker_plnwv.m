function tof_map = threshold_picker_plnwv( rfdata,kgrid )
%しきい値法で波形先頭を抽出する．
%初期化
dt = kgrid.dt;
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