function [ medium ] = makePatternMedium( pattern, param, rateEMCL, nEMCL, rateIMCL )
%MAKERANDOMMEDIUM： 割合，サイズを指定して複数条件のMEDIUMを生成する関数
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

%patternの決定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1:縦縞-高さ＝アレイサイズ（2018/08/05修正：rEMCLが不適切だったのを直した．アレイサイズに限定されていた
%                                              ことを考慮していなかった）
% 2:縦縞-高さ＝計算領域高さ
% 3:均質媒質
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pattern == 1
    center = zeros(1,nEMCL);
    unit.Nx = int16((Lx/ROI.x)*ROI.Nx*(rateEMCL/100)/nEMCL);%EMCLクラスタ一つにつき幅をどれだけにするか．
    offset_N = int16((Nx-ROI.Nx)/2);
    for i = 1:nEMCL
        center(1,i)= int16( i * (ROI.Nx / (nEMCL+1)));
        medium.density(offset_N+1:offset_N+ROI.Ny, offset_N+center(1,i)-unit.Nx/2+1:offset_N+center(1,i)+unit.Nx/2) = den_fat;
        medium.sound_speed(offset_N+1:offset_N+ROI.Ny, offset_N+center(1,i)-unit.Nx/2+1:offset_N+center(1,i)+unit.Nx/2) = v_fat;
    end
end

if pattern == 2
    center = zeros(1,nEMCL);
    unit.Nx = int16(ROI.Nx*(rateEMCL/100)/nEMCL);
    offset_N = int16((Nx-ROI.Nx)/2);
    for i = 1:nEMCL
        center(1,i)= int16( i * (ROI.Nx / (nEMCL+1)));
        medium.density(:, offset_N+center(1,i)-unit.Nx/2+1:offset_N+center(1,i)+unit.Nx/2) = den_fat;
        medium.sound_speed(:, offset_N+center(1,i)-unit.Nx/2+1:offset_N+center(1,i)+unit.Nx/2) = v_fat;
    end
end

end
end