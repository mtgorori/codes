% clear
% close all
% load param_ring.mat
%%%% ϊp[^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ή¬
v_water = 1540;%ΜΉ¬[m/s]
v_fat = 1420;%bΜΉ¬[m/s]
%§x
den_water = 1000;
den_fat = 920;
%bΜζ`σ
fat_cx = 0; fat_cy = 20.e-3; % SΚu
fat_radius = 10.e-3;
fd = fat_radius;
%Zέθ
Nx=int16(param.grid.Nx);
Ny=int16(param.grid.Ny);
kgrid = kWaveGrid(param.grid.Nx, param.grid.dx, param.grid.Ny, param.grid.dy);
x_vec=kgrid.x_vec;
y_vec=kgrid.y_vec;
%%%% Ή¬ͺz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
cd('\\Azlab-fs01\€Ί\Βlwork\|ΰ(Π)\USCTSim-master')
save('medium_03_cy_20','medium')