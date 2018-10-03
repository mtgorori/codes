function [ medium ] = makeAngleMedium( param, angle, rateEMCL, rateIMCL )
%�v�Z�̈�̒����ɑ΂��ā�2�{�̒����̗̈��ݒ�D�p�O�x�̔}�����쐬�D
%�摜�Ƃ��Ĉ����C��]�������|���C�Ō�Ɍv�Z�̈敪�g���~���O���������Ĕ}�������D
% ����
v_muscle = 1580;%�ؓ��̉���[m/s]
v_fat = 1450;%���b�̉���[m/s]
% ���x
den_muscle = 1040;%�ؓ��̖��x[kg/m3]
den_fat = 920;%���b�̖��x[kg/m3]
%medium�̏�����
Nx = param.grid.Nx;
Ny = param.grid.Ny;
dx = param.grid.dx;
dy = param.grid.dy;
Lx = Nx * dx;%�v�Z�̈��x�����̒���
Ly = Ny * dy;%�v�Z�̈��y�����̒���
ROI.x = 0.04;%ROI��x�����̒���
ROI.y = 0.04;%ROI��y�����̒���
ROI.Nx = int16(Nx * (ROI.x / Lx));%ROI��x�����̃O���b�h��
ROI.Ny = int16(Ny * (ROI.y / Ly));%ROI��y�����̃O���b�h��
Nx2 = ceil(Nx * sqrt(2));%�g���~���O�̂��߂ɗ̈�g��
Ny2 = ceil(Ny * sqrt(2)); %�g���~���O�̂��߂ɗ̈�g��
Lx2 = sqrt(2) * Lx;%�g���~���O�΍�
Ly2 = sqrt(2) * Ly;%�g���~���O�΍�
baseDensity = den_muscle + (den_fat - den_muscle)*rateIMCL/100;
baseSoundSpeed = v_muscle + (v_fat - v_muscle)*rateIMCL/100;
kgrid = kWaveGrid(Nx2, param.grid.dx, Ny2, param.grid.dy);
x_vec = kgrid.x_vec;
y_vec = kgrid.y_vec;
%medium���̔��f
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

