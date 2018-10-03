function [ medium ] = makeAngleMedium( param, angle, rateEMCL, rateIMCL )
%計算領域の長さに対して√2倍の長さの領域を設定．仰角０度の媒質を作成．
%画像として扱い，回転処理を掛け，最後に計算領域分トリミング処理をして媒質を作る．
% 音速
v_muscle = 1580;%筋肉の音速[m/s]
v_fat = 1450;%脂肪の音速[m/s]
% 密度
den_muscle = 1040;%筋肉の密度[kg/m3]
den_fat = 920;%脂肪の密度[kg/m3]
%mediumの初期化
Nx = param.grid.Nx;
Ny = param.grid.Ny;
dx = param.grid.dx;
dy = param.grid.dy;
Lx = Nx * dx;%計算領域のx方向の長さ
Ly = Ny * dy;%計算領域のy方向の長さ
ROI.x = 0.04;%ROIのx方向の長さ
ROI.y = 0.04;%ROIのy方向の長さ
ROI.Nx = int16(Nx * (ROI.x / Lx));%ROIのx方向のグリッド数
ROI.Ny = int16(Ny * (ROI.y / Ly));%ROIのy方向のグリッド数
Nx2 = ceil(Nx * sqrt(2));%トリミングのために領域拡大
Ny2 = ceil(Ny * sqrt(2)); %トリミングのために領域拡大
Lx2 = sqrt(2) * Lx;%トリミング対策
Ly2 = sqrt(2) * Ly;%トリミング対策
baseDensity = den_muscle + (den_fat - den_muscle)*rateIMCL/100;
baseSoundSpeed = v_muscle + (v_fat - v_muscle)*rateIMCL/100;
kgrid = kWaveGrid(Nx2, param.grid.dx, Ny2, param.grid.dy);
x_vec = kgrid.x_vec;
y_vec = kgrid.y_vec;
%medium情報の反映
% ROI.l = (ROI.y-sqrt(ROI.y^2-4*rateEMCL*ROI.y^2))/2;
ROI.l = (5*ROI.y-sqrt(25*ROI.y^2-20*rateEMCL*ROI.y^2))/2;
medium.density = ones(Ny2, Nx2) * baseDensity;
medium.sound_speed = ones(Ny2, Nx2) * baseSoundSpeed;
% center.Nx = Nx2;
% center.Ny = Ny2;
for i = 1:Nx2
    for j = 1:Ny2
        if (y_vec(i)<=ROI.l/2) && (-ROI.l/2<=y_vec(i)) && (x_vec(j)<=ROI.x/2-ROI.l/10)  && (-(ROI.x/2-ROI.l/10)<=x_vec(j))
            medium.sound_speed(j,i) = v_fat;
            medium.density(j,i) = den_fat;
        end
    end
end

 medium.sound_speed  = imrotate(medium.sound_speed,angle);
 [Ix, Iy] = size(medium.sound_speed);
 medium.sound_speed = medium.sound_speed(Ix/2-Nx/2:Ix/2+Nx/2-1,Iy/2-Ny/2:Iy/2+Ny/2-1);
 medium.sound_speed((~medium.sound_speed)) = baseSoundSpeed; 
 medium.density  = imrotate(medium.density,angle);
 [Ix2, Iy2] = size(medium.density);
 medium.density = medium.density(Ix2/2-Nx/2:Ix2/2+Nx/2-1,Iy2/2-Ny/2:Iy2/2+Ny/2-1);
 medium.density((~medium.density)) = baseDensity; 
end

