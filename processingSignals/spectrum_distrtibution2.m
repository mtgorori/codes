function spectrum_distrtibution2( case_number,medium_path,save_path )
close all

cd('H:data/configration')
load('transducer_2board.mat')

%�����ݒ�
cd(medium_path);
load('case1.mat','kgrid')
x_vec = kgrid.x_vec(51:450);
y_vec = kgrid.y_vec(51:450);
    
    interp_x = 1:5000;
    interp_y = 1:5000;
x_vec = interp1(x_vec,interp_x);
leng_x = 5000;
leng_y = 5000;
theta = zeros(5000,5000);
rho = zeros(5000,5000);
interval_num = 360;
reference_angle = linspace(0, ((interval_num-1)/interval_num)*2*pi, interval_num)';
d_theta = pi/interval_num;
spectrum_distribution = zeros(interval_num,case_number);
FFT = zeros(5000,5000);
scan_domain = zeros(5000,5000);
target_rho = (rho<=min(max(x_vec)*1000,max(y_vec)*1000));
major_angle = zeros(1,case_number);
for i = 1:leng_y
    for j = 1:leng_x
        [theta(i,j),rho(i,j)] = cart2pol(x_vec(j)*1000,y_vec(i)*1000);
        if theta(i,j) < 0
            theta(i,j) = theta(i,j) + 2*pi;
        end
    end
end


%�f�[�^���
for ii = 1:case_number
    cd(medium_path);
    loadfilename =  sprintf('case%d.mat',ii);
    load(loadfilename);
    cd(save_path);
    dirname = sprintf('case%d',ii);
    mkdir(dirname);
    
    figure(1);
    imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
    xlabel('x��[mm]')
    ylabel('y��[mm]')
    colorbar;
    caxis([1450 1580]);
    c = colorbar;
    c.Label.String = '[m/s]';
    ax = gca;
    ax.YDir = 'normal';
    axis equal
    axis tight
    myfilename0 = sprintf('/case%d/case%d_SSdist',ii,ii);
    myfilename = strcat(save_path,myfilename0);
    exportfig(myfilename,'png',[300,300]);
    
    figure(2);
    tmp = medium.sound_speed(51:450,51:450);
    binary_medium = (tmp == min(min(tmp)));
    imagesc(kgrid.x_vec(51:450)*1000,kgrid.y_vec(51:450)*1000,binary_medium);
    colorbar;
    axis equal
    axis tight
    xlabel('x��[mm]')
    ylabel('y��[mm]')
    ax = gca;
    ax.YDir = 'normal';
    myfilename0 = sprintf('/case%d/case%d_binalize',ii,ii);
    myfilename = strcat(save_path,myfilename0);
    exportfig(myfilename,'png',[300,300]);
    
    figure(3);
    fft2_binary_medium = fft2(binary_medium);
    imagesc(log(abs(fftshift(fft2_binary_medium))));
    colorbar;
    ax = gca;
    ax.YDir = 'normal';
    xlabel('�������g��')
    ylabel('�������g��')
    axis equal
    axis tight
    axis off
    myfilename0 = sprintf('/case%d/case%d_fft_log',ii,ii);
    myfilename = strcat(save_path,myfilename0);
    exportfig(myfilename,'png',[300,300]);
    
    interp_x = 1:5000;
    interp_y = 1:5000;
    FFT(:,:)  = interp2(abs(fftshift(fft2_binary_medium)),interp_y,interp_x);
    
    target_theta = ((2*pi-d_theta<=theta)&(theta<=2*pi))|((0<=theta)&(theta<d_theta));
    scan_domain(:,:) = target_rho.*target_theta;
    spectrum_distribution(1,ii) = sum(sum(scan_domain.*FFT(:,:)))/(sum(sum(scan_domain)));
    
    
    for jj = 2:interval_num
        target_theta = (reference_angle(jj,1)-d_theta<theta)&(theta<=reference_angle(jj,1)+d_theta);
        scan_domain = target_rho.*target_theta;
        spectrum_distribution(jj,ii) = sum(sum(scan_domain.*FFT(:,:)))/(sum(sum(scan_domain)));
    end
    
    figure(4);
    plot(reference_angle,spectrum_distribution(:,ii));
    xlabel('angle[rad]')
    ylabel('power')
    xlim([-0.5,2*pi+0.5])
    myfilename0 = sprintf('/case%d/case%d_spectrum_distribution',ii,ii);
    myfilename = strcat(save_path,myfilename0);
    exportfig(myfilename,'png',[300,300]);
    
    
    %%movie
    clear mov
    fr(1:interval_num) = struct('cdata',[],'colormap',[]);
    figure(5);
    % make movie
    for i=1:interval_num
        subplot(2,2,1);
        imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,binary_medium);
        colorbar;
        axis equal
        axis tight
        xlabel('x��[mm]')
        ylabel('y��[mm]')
        ax = gca;
        ax.YDir = 'normal';
        subplot(2,2,2);
        imagesc(20*log10(abs(FFT(:,:))));
        axis equal
        axis tight
        ax = gca;
        ax.YDir = 'normal';
        hold on
        image(100*scan_domain(:,:,i),'AlphaData',0.4);
        colorbar;
        subplot(2,2,[3,4]);
        plot(reference_angle(1:i),spectrum_distribution(1:i,ii));
        hold on
        scatter(reference_angle(i),spectrum_distribution(i,ii));
        xlabel('angle[rad]')
        ylabel('power')
        xlim([-0.5,2*pi+0.5])
        pause(0.05);
        hold off
        drawnow;
        fr(i) = getframe(5);
    end
    % play movie
%     figure;
%     movie(fr)
    % save movie
    myfilename0 = sprintf('/case%d/case%d_spectrum_movie',ii,ii);
    myfilename = strcat(save_path,myfilename0);
    mv = VideoWriter(myfilename,'MPEG-4'); %#ok<TNMLP>
    mv.FrameRate = 36; % �� fps�Ɠ��� %ART:3
    open(mv);
    writeVideo(mv,fr);
    close(mv);
    
    %������p�x����
    [~,ind_spectrum] = max(spectrum_distribution(:,ii));
    major_angle(1,ii) = reference_angle(ind_spectrum);
    
end

myfilename0 = sprintf('work_space');
myfilename = strcat(save_path,myfilename0);
save(myfilename,'major_angle','reference_angle','spectrum_distribution')

end

