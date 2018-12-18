% clear
% close all
% load param_ring.mat
%%%% �����p�����[�^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
v_water = 1540;%���̉���[m/s]
v_fat = 1420;%���b�̉���[m/s]
%���x
den_water = 1000;
den_fat = 920;
%���b�̈�`��
fat_cx = 0; fat_cy = 0; % ���S�ʒu
fat_radius = 10.e-3;
fd = fat_radius;
%�Z���ݒ�
Nx=int16(param.grid.Nx);
Ny=int16(param.grid.Ny);
kgrid = kWaveGrid(param.grid.Nx, param.grid.dx, param.grid.Ny, param.grid.dy);
x_vec=kgrid.x_vec;
y_vec=kgrid.y_vec;
%%%% �������z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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