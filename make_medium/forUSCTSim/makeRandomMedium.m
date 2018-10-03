function [ medium ] = makeRandomMedium( param, rateEMCL, nEMCL, rateIMCL )
%MAKERANDOMMEDIUM�F �����C�T�C�Y���w�肵�ĕ���������MEDIUM�𐶐�����֐�
% ����
v_muscle = 1580;%�ؓ��̉���[m/s]
v_fat = 1450;%���b�̉���[m/s]
% ���x
den_muscle = 1040;%�ؓ��̖��x[kg/m3]
den_fat = 920;%���b�̖��x[kg/m3]
%medium�̏�����
Nx= param.grid.Nx;
Ny= param.grid.Ny;
dx = param.grid.dx;
dy = param.grid.dy;
baseDensity = den_muscle + (den_fat - den_muscle)*rateIMCL/100;
baseSoundSpeed = v_muscle + (v_fat - v_muscle)*rateIMCL/100;
medium.density = ones(Ny, Nx) * baseDensity;
medium.sound_speed = ones(Ny, Nx)* baseSoundSpeed;
%kgrid�쐬
kgrid = kWaveGrid(param.grid.Nx, param.grid.dx, param.grid.Ny, param.grid.dy);
x_vec=kgrid.x_vec;
y_vec=kgrid.y_vec;
%EMCL�̃h���C���ݒ�i���a������o���j
totalArea = (Nx*dx) * (Ny*dy);
EMCLArea = (rateEMCL/100) * totalArea;
radEMCL2 = (EMCLArea/nEMCL)/pi;
%EMCL�̒��S���W�z��̐ݒ�
numCircles = 1;
numCirclesMax = nEMCL+1;
maxIterations = 500 * numCirclesMax; % Fail Safe.
iteration = 1;
% r = 1000000;
% while abs(r*100 - rateEMCL)  > 0.3
    while  numCircles < numCirclesMax && iteration < maxIterations
        xTrial = x_vec(randi(Nx));
        yTrial = y_vec(randi(Ny));
        iteration = iteration + 1; % Fail safe.
        % Only need to check for overlap for second and later circles.
        if numCircles > 1
            % Find distance from other, prior circles.
            distances = sqrt((xTrial - x) .^ 2 + (yTrial - y) .^ 2);
            if min(distances) < 2 * sqrt(radEMCL2)
                % It's overlapping at least one of the prior ones
                continue; % Skip to end of loop and continue with loop.
                %if
            end
        end
        x(numCircles) = xTrial;
        y(numCircles) = yTrial;
        numCircles = numCircles + 1;
    end
    radii = radEMCL2 * ones(1, length(x));
    for n = 1:length(radii)
        for i = 1:Nx
            for j = 1:Ny
                if ((x_vec(i)-x(n))^2 + (y_vec(j)-y(n))^2)<= radii(n)%radii(n)^2
                    medium.sound_speed(j,i) = v_fat;
                    medium.density(j,i) = den_fat;
                end
            end
        end
    end
%     r = nnz(medium.sound_speed - v_fat);
%     r = 1-r/(Nx*Ny);
% end
end

