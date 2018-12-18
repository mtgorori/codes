[us,RMSE] = SART_reconstruction(projection);
cd('../result');
save('SART_aic_window150','us','RMSE');
clear 
cd('../USCTSim-master/result/kwave_2board_object');
load('rfdata.mat');
rfdata_ob = rfdata;
cd('../kwave_2board_water');
load('rfdata.mat');
load('kgrid.mat');
load('sensor.mat');
rfdata_wa = rfdata;
clear rfdata;
cd('H:codes');
tof_map = get_tof_AIC(rfdata_ob,rfdata_wa,kgrid.mat,sensor.mat,75);
cd('H:result');  
save('tof_map_win_150_right_shift_2board.mat','tof_map');
projection = tof_map;
cd('H:codes');
[us,RMSE] = SART_reconstruction_2board(projection);
cd('H:result'); 
save('SART_recon_2board_right_shifted_20mm_kwave.mat','us','RMSE');
clear


