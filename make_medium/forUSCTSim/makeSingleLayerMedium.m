function [ medium ] = makeSingleLayerMedium( param, kgrid, rateIMCL, x_coordinate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[ medium ] = makeSingleLayerMedium( param, kgrid, rateIMCL, x_coordinate)
%点媒質を作成する．
%x_coordinate:[x1,x2,...,xn](mm)
%y_coordinate:[y1,y2,...,yn](mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 初期パラメタ設定
% 音速
v_muscle = 1580;%筋肉の音速[m/s]
v_fat = 1450;%脂肪の音速[m/s]
% 密度
den_muscle = 1040;%筋肉の密度[kg/m3]
den_fat = 920;%脂肪の密度[kg/m3]
% alpha_coeff_muscle = 0.74;
% alpha_coeff_fat = 0.29;
% BonA_muscle = 6.6;
% BonA_fat = 10.0;
%mediumの初期化
Nx = param.grid.Nx;
Ny = param.grid.Ny;
dx = param.grid.dx;
dy = param.grid.dy;
Lx = Nx * dx;%全領域のx方向の長さ
Ly = Ny * dy;%全領域のy方向の長さ
ROI.x = 0.04;%ROIのx方向の長さ
ROI.y = 0.04;%ROIのy方向の長さ
ROI.Nx = int16(Nx * (ROI.x / Lx));%ROIのx方向のグリッド数
ROI.Ny = int16(Ny * (ROI.y / Ly));%ROIのy方向のグリッド数
baseDensity = den_muscle + (den_fat - den_muscle)*rateIMCL/100;
baseSoundSpeed = v_muscle + (v_fat - v_muscle)*rateIMCL/100;
medium.density = ones(Ny, Nx) * baseDensity;
medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
% 層媒質設定
tmp_x = find(kgrid.x_vec == x_coordinate/1e3);
medium.sound_speed(:,1:tmp_x) = v_fat;
medium.density(:,1:tmp_x) = den_fat;
end

