function [ medium ] = makeSingleLayerMedium( param, kgrid, rateIMCL, x_coordinate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[ medium ] = makeSingleLayerMedium( param, kgrid, rateIMCL, x_coordinate)
%�_�}�����쐬����D
%x_coordinate:[x1,x2,...,xn](mm)
%y_coordinate:[y1,y2,...,yn](mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% �����p�����^�ݒ�
% ����
v_muscle = 1580;%�ؓ��̉���[m/s]
v_fat = 1450;%���b�̉���[m/s]
% ���x
den_muscle = 1040;%�ؓ��̖��x[kg/m3]
den_fat = 920;%���b�̖��x[kg/m3]
% alpha_coeff_muscle = 0.74;
% alpha_coeff_fat = 0.29;
% BonA_muscle = 6.6;
% BonA_fat = 10.0;
%medium�̏�����
Nx = param.grid.Nx;
Ny = param.grid.Ny;
dx = param.grid.dx;
dy = param.grid.dy;
Lx = Nx * dx;%�S�̈��x�����̒���
Ly = Ny * dy;%�S�̈��y�����̒���
ROI.x = 0.04;%ROI��x�����̒���
ROI.y = 0.04;%ROI��y�����̒���
ROI.Nx = int16(Nx * (ROI.x / Lx));%ROI��x�����̃O���b�h��
ROI.Ny = int16(Ny * (ROI.y / Ly));%ROI��y�����̃O���b�h��
baseDensity = den_muscle + (den_fat - den_muscle)*rateIMCL/100;
baseSoundSpeed = v_muscle + (v_fat - v_muscle)*rateIMCL/100;
medium.density = ones(Ny, Nx) * baseDensity;
medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
% �w�}���ݒ�
tmp_x = find(kgrid.x_vec == x_coordinate/1e3);
medium.sound_speed(:,1:tmp_x) = v_fat;
medium.density(:,1:tmp_x) = den_fat;
end

