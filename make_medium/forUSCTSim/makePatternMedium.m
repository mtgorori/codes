function [ medium ] = makePatternMedium( pattern, param, rateEMCL, nEMCL, rateIMCL )
%MAKERANDOMMEDIUM�F �����C�T�C�Y���w�肵�ĕ���������MEDIUM�𐶐�����֐�
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

%pattern�̌���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1:�c��-�������A���C�T�C�Y�i2018/08/05�C���FrEMCL���s�K�؂������̂𒼂����D�A���C�T�C�Y�Ɍ��肳��Ă���
%                                              ���Ƃ��l�����Ă��Ȃ������j
% 2:�c��-�������v�Z�̈捂��
% 3:�ώ��}��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pattern == 1
    center = zeros(1,nEMCL);
    unit.Nx = int16((Lx/ROI.x)*ROI.Nx*(rateEMCL/100)/nEMCL);%EMCL�N���X�^��ɂ������ǂꂾ���ɂ��邩�D
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