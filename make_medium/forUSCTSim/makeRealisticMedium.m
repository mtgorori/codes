function [ medium, rateEMCL ] = makeRealisticMedium( param, rateIMCL, num_pict,location,max_num_pict,lossless)
%% �����p�����^�ݒ�
% ����
v_muscle = 1580;%�ؓ��̉���[m/s]
v_fat = 1450;%���b�̉���[m/s]

% ���x
den_muscle = 1040;%�ؓ��̖��x[kg/m3]
den_fat = 920;%���b�̖��x[kg/m3]
alpha_coeff_muscle = 0.74;
alpha_coeff_fat = 0.29;
BonA_muscle = 6.6;
BonA_fat = 10.0;
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

%���������l�������Ŏw��(EMCL��L���̃����W24~36 %���x)
if (1<= location) &&  (location<= 3)
    cd('H:\data\pictures_of_calf\main_data\2018_07_10');
    I = imread('IMG_2918.jpg');
elseif (4 <= location) && (location <= 8)
    cd('H:\data\pictures_of_calf\main_data\2018_08_28');
    myfilename = sprintf('%d.jpg',location);
    I = imread(myfilename);
    %���������l���}�j���A���Ŏw��(EMCL��L���̃����W10~22 %���x)
elseif (9 <= location ) && (location <= 11)
    cd('H:\data\pictures_of_calf\main_data\2018_07_10');
    I = imread('IMG_2918.jpg');
elseif (12 <= location) && (location <= 16)
    cd('H:\data\pictures_of_calf\main_data\2018_08_28');
    myfilename = sprintf('%d.jpg',location-8);
    I = imread(myfilename);
end

if location == 1
    frame_x15 = 3600;%x�����̃t���[����(15 cm����)
    frame_x5 = frame_x15/3;%5cm����=1200
    frame_y4 = frame_x5 * 4/5;
    min_x = 2000;
    min_y = 1300;
    max_x = 3200;
    max_y = 3000;
    [I2, ~] = imcrop(I,[min_x (min_y+(num_pict-1)*(max_y - min_y - frame_y4)/(max_num_pict-1)) frame_x5-1 frame_y4-1]);
    I2 = imrotate(I2,90);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 2
    frame_x15 = 3600;%x�����̃t���[����(15 cm����)
    frame_x5 = frame_x15/3;%5cm����=1200
    frame_y4 = frame_x5 * 4/5;
    min_x = 200;
    min_y = 1100;
    max_x = 1600;
    max_y = 3200;
    [I2, ~] = imcrop(I,[min_x (min_y+(num_pict-1)*(max_y - min_y - frame_y4)/(max_num_pict-1)) frame_x5-1 frame_y4-1]);
    I2 = imrotate(I2,90);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 3
    frame_y15 = 3600;%x�����̃t���[����(15 cm����)
    frame_y5 = frame_y15/3;%5cm����=1200
    frame_x4 = frame_y5 * 4/5;
    min_x = 2250;
    min_y = 1800;
    max_x = 4800;
    max_y = 3000;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 4
    frame_y5 = 650;
    frame_x4 = frame_y5 * 4/5;
    min_x = 245;
    min_y = 820;
    max_x = 2420;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end
if location == 5
    frame_y5 = 840;
    frame_x4 = frame_y5 * 4/5;
    min_x = 1185;
    min_y = 547;
    max_x = 3340;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end
if location == 6
    frame_y5 = 800;
    frame_x4 = frame_y5 * 4/5;
    min_x = 213;
    min_y = 595;
    max_x = 2385;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end
if location == 7
    frame_y5 = 850;
    frame_x4 = frame_y5 * 4/5;
    min_x = 1270;
    min_y = 678;
    max_x = 3400;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 8
    frame_y5 = 900;
    frame_x4 = frame_y5 * 4/5;
    min_x = 378;
    min_y = 525;
    max_x = 2475;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end
if location == 9
    frame_x15 = 3600;%x�����̃t���[����(15 cm����)
    frame_x5 = frame_x15/3;%5cm����=1200
    frame_y4 = frame_x5 * 4/5;
    min_x = 2000;
    min_y = 1300;
    max_x = 3200;
    max_y = 3000;
    [I2, ~] = imcrop(I,[min_x (min_y+(num_pict-1)*(max_y - min_y - frame_y4)/(max_num_pict-1)) frame_x5-1 frame_y4-1]);
    I2 = imrotate(I2,90);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 10
    frame_x15 = 3600;%x�����̃t���[����(15 cm����)
    frame_x5 = frame_x15/3;%5cm����=1200
    frame_y4 = frame_x5 * 4/5;
    min_x = 200;
    min_y = 1100;
    max_x = 1600;
    max_y = 3200;
    [I2, ~] = imcrop(I,[min_x (min_y+(num_pict-1)*(max_y - min_y - frame_y4)/(max_num_pict-1)) frame_x5-1 frame_y4-1]);
    I2 = imrotate(I2,90);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 11
    frame_y15 = 3600;%x�����̃t���[����(15 cm����)
    frame_y5 = frame_y15/3;%5cm����=1200
    frame_x4 = frame_y5 * 4/5;
    min_x = 2250;
    min_y = 1800;
    max_x = 4800;
    max_y = 3000;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end
if location == 12
    frame_y5 = 650;
    frame_x4 = frame_y5 * 4/5;
    min_x = 245;
    min_y = 820;
    max_x = 2420;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh*1.15);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
%     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
%     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
%     medium.density = ones(Ny, Nx) * baseDensity;
%     medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
%     baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
%     baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
%     medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
%     medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 13
    frame_y5 = 840;
    frame_x4 = frame_y5 * 4/5;
    min_x = 1185;
    min_y = 547;
    max_x = 3340;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh*1.12);%�������l����2018/10/12
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 14
    frame_y5 = 800;
    frame_x4 = frame_y5 * 4/5;
    min_x = 213;
    min_y = 595;
    max_x = 2385;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh*1.08);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 15
    frame_y5 = 850;
    frame_x4 = frame_y5 * 4/5;
    min_x = 1270;
    min_y = 678;
    max_x = 3400;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh*1.17);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.1;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

if location == 16
    frame_y5 = 900;
    frame_x4 = frame_y5 * 4/5;
    min_x = 378;
    min_y = 525;
    max_x = 2475;
    [I2, ~] = imcrop(I,[(min_x+(num_pict-1)*(max_x - min_x - frame_x4)/(max_num_pict-1)) min_y frame_x4-1 frame_y5-1]);
    grayImg = rgb2gray(I2);
    thresh = graythresh(grayImg);
    bwImg_pre = imbinarize(grayImg, thresh*1.13);
    bwImg = imresize(bwImg_pre, [Ny ROI.Nx]);
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
    if lossless == 0
        medium.BonA = ones(Ny,Nx) * BonA_muscle;
        medium.alpha_coeff = ones(Ny,Nx) * alpha_coeff_muscle;
        medium.alpha_power = 1.05;
        medium.alpha_coeff(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * alpha_coeff_fat + ~bwImg * alpha_coeff_muscle;
        medium.BonA(:,Nx/2-ROI.Nx/2:Nx/2+ROI.Nx/2) = bwImg * BonA_fat + ~bwImg * BonA_muscle;
        %�v�C��2018/09/28
    end
    %     rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) /
    %     numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2,:))) * 100;%2018/10/09�ȑO
    rateEMCL = (sum(sum(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) / numel(bwImg(Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2+1,:))) * 100;
    medium.density = ones(Ny, Nx) * baseDensity;
    medium.sound_speed = ones(Ny, Nx) * baseSoundSpeed;
    baseDensity = (den_muscle + (den_fat - den_muscle)*rateIMCL/(100-rateEMCL));
    baseSoundSpeed = (v_muscle + (v_fat - v_muscle)*rateIMCL/(100-rateEMCL));
    medium.density(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * den_fat + ~bwImg * baseDensity;
    medium.sound_speed(:,Nx/2-ROI.Nx/2+1:Nx/2+ROI.Nx/2) = bwImg * v_fat + ~bwImg * baseSoundSpeed;
end

end