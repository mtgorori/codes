% clear
% close all
% load param_ring.mat
%%%% 初期パラメータ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%音速
v_water = 1540;%水の音速[m/s]
v_fat = 1420;%脂肪の音速[m/s]
%密度
den_water = 1000;
den_fat = 920;
%脂肪領域形状
fat_cx = 0; fat_cy = 0; % 中心位置
fat_radius = 10.e-3;
fd = fat_radius;
%セル設定
Nx=int16(param.grid.Nx);
Ny=int16(param.grid.Ny);
kgrid = kWaveGrid(param.grid.Nx, param.grid.dx, param.grid.Ny, param.grid.dy);
x_vec=kgrid.x_vec;
y_vec=kgrid.y_vec;
%%%% 音速分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
medium.sound_speed = v_water*ones(Nx,Ny);
medium.density=den_water*ones(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if ((x_vec(i)-fat_cx)^2 + (y_vec(j)-fat_cy)^2)<= fd^2
            medium.sound_speed(i,j) = v_fat;
            medium.density(i,j) = den_fat;
        end
    end
end
figure;
imagesc(medium.sound_speed);
figure;
imagesc(medium.density);
cd('H:USCTSim-master')
save('medium_04_center','medium')