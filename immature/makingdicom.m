% clear
region_nx = 1024;
region_ny = region_nx;
path1='G:\GdriveBackup0707\kondo\SEKI MIKA';
path2='1.3.46.670589.11.71301.5.0.600.2015120713070416000';
% path2='1.3.46.670589.11.71301.5.0.9424.2015120713065106000';

searching = 1;  % 1だとフォルダ内検索、0だと単体
%     fname='1.3.46.670589.11.71301.5.0.9424.2015120713065160004';
    fnames=dir(strcat(path1,'\',path2,'\*.dcm'));
    data=zeros(1024,1024,length(fnames));
if searching ==1
for fnum=1:length(fnames)
mridata_dicom = dicomread(strcat(path1,'\',path2,'\',fnames(fnum).name));
[X,Y] = meshgrid(1:512);

% 1024メッシュの場合
[Xq,Yq] = meshgrid(0.5:0.5:512);
mridata_dicom1024 = interp2(X,Y,single(mridata_dicom),Xq,Yq);
adipose_region = mridata_dicom1024>1000;
water_region = mridata_dicom1024<100;
skin_region = mridata_dicom1024<450&mridata_dicom1024>100;
skin_region(394:620,430:618) = 0;
skin_region(628:666,486:602) = 0;
skin_region(526,408:409) = 0;
glandular_region = ones(1024) - skin_region - water_region - adipose_region;

% 2048メッシュの場合
% [Xq,Yq] = meshgrid(0.25:0.25:512);
% mridata_dicom2048 = interp2(X,Y,single(mridata_dicom),Xq,Yq);
% adipose_region = mridata_dicom2048>1000;
% water_region = mridata_dicom2048<100;
% skin_region = mridata_dicom2048<450&mridata_dicom2048>100;
% skin_region(394*2:620*2,430*2:618*2) = 0;
% skin_region(628*2:666*2,486*2:602*2) = 0;
% skin_region(526*2,408*2:409*2) = 0;
% glandular_region = ones(2048) - skin_region - water_region - adipose_region;

% 4096メッシュの場合
% [Xq,Yq] = meshgrid(0.125:0.125:512);
% mridata_dicom4096 = interp2(X,Y,single(mridata_dicom),Xq,Yq);
% adipose_region = mridata_dicom4096>1000;
% water_region = mridata_dicom4096<100;
% skin_region = mridata_dicom4096<450&mridata_dicom4096>100;
% skin_region(394*4:620*4,430*4:618*4) = 0;
% skin_region(628*4:666*4,486*4:602*4) = 0;
% skin_region(526*4,408*4:409*4) = 0;
% glandular_region = ones(4096) - skin_region - water_region - adipose_region;


wire_region = zeros(region_nx,region_ny,'logical');
% wire_region(region_nx/2,region_ny/2) = 1;
% xx_wire_num = 6;
% yy_wire_num = 1;
% for xxwi = 1:xx_wire_num
%     for yywi = 1:yy_wire_num
% % %         glandular_region = glandular_region + makeDisc(Nx,Ny,offset+xxgl*Nx/16,offset+yygl*Ny/16,8);
%          wire_region(region_nx/2 + (xxwi-1)*region_nx/32, region_ny/2) = 1;
%          glandular_region(region_nx/2 + (xxwi-1)*region_nx/32, region_ny/2) = 0;
%          adipose_region(region_nx/2 + (xxwi-1)*region_nx/32, region_ny/2) = 0;
%     end
% end
glandular_region = logical(glandular_region);

% figure;imagesc(skin_region+glandular_region*2+adipose_region*3+wire_region*4);title(fnames(fnum).name)
data(:,:,fnum)=skin_region+glandular_region*2+adipose_region*3+wire_region*4;

disp(['Finished num: ',num2str(fnum),'/',num2str(length(fnames))])

end
else
    mridata_dicom = dicomread(strcat(path1,'\',path2,'\',fname));
[X,Y] = meshgrid(1:512);

% 1024メッシュの場合
[Xq,Yq] = meshgrid(0.5:0.5:512);
mridata_dicom1024 = interp2(X,Y,single(mridata_dicom),Xq,Yq);
adipose_region = mridata_dicom1024>1000;
water_region = mridata_dicom1024<100;
skin_region = mridata_dicom1024<450&mridata_dicom1024>100;
skin_region(394:620,430:618) = 0;
skin_region(628:666,486:602) = 0;
skin_region(526,408:409) = 0;
glandular_region = ones(1024) - skin_region - water_region - adipose_region;

% 2048メッシュの場合
% [Xq,Yq] = meshgrid(0.25:0.25:512);
% mridata_dicom2048 = interp2(X,Y,single(mridata_dicom),Xq,Yq);
% adipose_region = mridata_dicom2048>1000;
% water_region = mridata_dicom2048<100;
% skin_region = mridata_dicom2048<450&mridata_dicom2048>100;
% skin_region(394*2:620*2,430*2:618*2) = 0;
% skin_region(628*2:666*2,486*2:602*2) = 0;
% skin_region(526*2,408*2:409*2) = 0;
% glandular_region = ones(2048) - skin_region - water_region - adipose_region;

% 4096メッシュの場合
% [Xq,Yq] = meshgrid(0.125:0.125:512);
% mridata_dicom4096 = interp2(X,Y,single(mridata_dicom),Xq,Yq);
% adipose_region = mridata_dicom4096>1000;
% water_region = mridata_dicom4096<100;
% skin_region = mridata_dicom4096<450&mridata_dicom4096>100;
% skin_region(394*4:620*4,430*4:618*4) = 0;
% skin_region(628*4:666*4,486*4:602*4) = 0;
% skin_region(526*4,408*4:409*4) = 0;
% glandular_region = ones(4096) - skin_region - water_region - adipose_region;


wire_region = zeros(region_nx,region_ny,'logical');
% wire_region(region_nx/2,region_ny/2) = 1;
% xx_wire_num = 6;
% yy_wire_num = 1;
% for xxwi = 1:xx_wire_num
%     for yywi = 1:yy_wire_num
% % %         glandular_region = glandular_region + makeDisc(Nx,Ny,offset+xxgl*Nx/16,offset+yygl*Ny/16,8);
%          wire_region(region_nx/2 + (xxwi-1)*region_nx/32, region_ny/2) = 1;
%          glandular_region(region_nx/2 + (xxwi-1)*region_nx/32, region_ny/2) = 0;
%          adipose_region(region_nx/2 + (xxwi-1)*region_nx/32, region_ny/2) = 0;
%     end
% end
glandular_region = logical(glandular_region);

% figure;imagesc(skin_region+glandular_region*2+adipose_region*3+wire_region*4);title(fnames(fnum).name)
data=skin_region+glandular_region*2+adipose_region*3+wire_region*4;
end
% clearvars -except data
% plotcheck
% save('region_map_mri1024.mat');