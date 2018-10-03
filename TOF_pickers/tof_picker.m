%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          tof_picker
%   Description:    several kinds of TOF pickers 
%   Author:         Xiaolei Qu
%   Date:           2014.06.13
%   Version:        1.0
%   History:        1.0  threshold(TH) and cross correlation(CC) methods are writtern basing on previous code.
%
%   Input:
%   When input parameter number is 3 (nargin == 3)
%
%       varargin{1}     receive_data:   All receive data during UCT simulation. 3D data, 
%                       1st D. is for time, 
%                       2nd D. is for receive element, 
%                       3rd D. is for transmit element.
%       varargin{2}     time_hist:      Receive data time mark(or time hist). 2D data, 
%                       1st D. is for time, 
%                       2nd D. is for transmit element. 
%                       All receive element share same time mark for same transmit element.
%       varargin{3}     method_type:    
%                       TH: threshold method
%                       AIC: Akaike information criterion method.
%                       WAIC: weighted AIC method.
%                       AICCC: combination of AIC and neighbor cross correlation.
%
%   When input parameter number is 5 (nargin == 5)
%
%       varargin{1}     receive_data for pattern:   All receive data during UCT simulation. 3D data, 
%                       1st D. is for time, 
%                       2nd D. is for receive element, 
%                       3rd D. is for transmit element.
%       varargin{2}     time_hist for pattern:      Receive data time mark(or time hist). 2D data, 
%                       1st D. is for time, 
%                       2nd D. is for transmit element. 
%                       All receive element share same time mark for same transmit element.
%       varargin{3}     similar to varargin{1} but receive data for water
%       varargin{4}     similar to varargin{2} but time hist for water
%       varargin{5}     method_type: (a string)    
%                       TH: threshold method
%                       AIC: Akaike information criterion method.
%                       WAIC: weighted AIC method.
%                       AICCC: combination of AIC and neighbor cross correlation.
%                       CC: cross correlation method (between water and pattern signal)
% 
%   Output:
%       tof_map:        time of flight map
%                       1st D. is for receive elements
%                       2nd D. is for transmit elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tof_map = tof_picker(varargin)         % main function for this file.
    
    %important parameters
    speed_in_water = 1491;  %
    ring_diameter = 0.2;    % 3.5cm diameter of Ring


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculated one single time of flight map
    % varargin{1} is the receive data
    % varargin{2} is the time hist of receive data
    % varargin{3} is the method typy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 3   % calculated one tof map
        
        % threshold method for TOF
        if strcmp(varargin{3},'TH')
            display('start running threshold method for TOF');
            tic;
            tof_map = tof_picker_TH(varargin{1}, varargin{2});
            display(['threshold method has gotten TOF and used ',num2str(toc),' s']);
            
        % first peak method
        elseif strcmp(varargin{3},'FP')
            display('start running FP method for TOF');
            tic;
            tof_map = tof_picker_FP(varargin{1}, varargin{2});
            display(['FP method has gotten TOF and used ',num2str(toc),' s']);
        
        % AIC method for TOF
        elseif strcmp(varargin{3},'AIC')
            display('start running AIC method for TOF');
            tic;
            tof_map = tof_picker_AIC(varargin{1}, varargin{2}, speed_in_water, ring_diameter);
            display(['AIC method has gotten TOF and used ',num2str(toc),' s']);
            
        % WAIC method for TOF
        elseif strcmp(varargin{3},'WAIC')
            display('start running WAIC method for TOF');
            tic;
            tof_map = tof_picker_WAIC(varargin{1}, varargin{2}, speed_in_water, ring_diameter);
            display(['WAIC method has gotten TOF and used ',num2str(toc),' s']);
        
        % AICC method for TOF
        % AICCC average method
        elseif strcmp(varargin{3},'AICNCC')
            display('start running AICNCC method for TOF');
            tic;
            tof_map = tof_picker_AICNCC(varargin{1}, varargin{2}, speed_in_water, ring_diameter);
            display(['AICNCC method has gotten TOF and used ',num2str(toc),' s']);
            
        else
            display(['there is no the method, ',varargin{3},', please checke the method name!!!']);
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculated tof difference which is (pattern tof - water tof)
    % varargin{1} is the receive data for pattern 
    % varargin{2} is the time hist for pattern 
    % varargin{3} is the receive data for water
    % varargin{4} is the time hist for water
    % varargin{3} is the method typy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    elseif nargin == 5  % calculated difference of two tof map
        
        % Threshold method for TOF difference
        if strcmp(varargin{5},'TH')
            display('start running threshold method for TOF difference');
            tic;
            tof_map = tof_picker_TH(varargin{1}, varargin{2}) - tof_picker_TH(varargin{3}, varargin{4});
            display(['threshold method has gotten TOF diff and used ',num2str(toc),' s']);
        
        % first peak method. 
        elseif strcmp(varargin{5},'FP')
            display('start running FP method for TOF difference');
            tic;
            tof_map = tof_picker_FP(varargin{1}, varargin{2}) - tof_picker_FP(varargin{3}, varargin{4});
            display(['FP method has gotten TOF diff and used ',num2str(toc),' s']);
            
        % AIC method for TOF difference
        elseif strcmp(varargin{5},'AIC')
            display('start running AIC method for TOF difference');
            tic;
            tof_map = tof_picker_AIC(varargin{1}, varargin{2}, speed_in_water, ring_diameter) - tof_picker_AIC(varargin{3}, varargin{4}, speed_in_water, ring_diameter);
            display(['AIC method has gotten TOF diff and used ',num2str(toc),' s']);
            
        % WAIC method for TOF difference
        elseif strcmp(varargin{5},'WAIC')
            display('start running WAIC method for TOF difference');
            tic;
            tof_map = tof_picker_WAIC(varargin{1}, varargin{2}, speed_in_water, ring_diameter) - tof_picker_WAIC(varargin{3}, varargin{4}, speed_in_water, ring_diameter);
            display(['WAIC method has gotten TOF diff and used ',num2str(toc),' s']);
            
        % AICNCC method for TOF difference
        elseif strcmp(varargin{5},'AICNCC')
            display('start running AICNCC method for TOF difference');
            tic;
            tof_map = tof_picker_AICNCC(varargin{1}, varargin{2}, speed_in_water, ring_diameter) - tof_picker_AICNCC(varargin{3}, varargin{4}, speed_in_water, ring_diameter);
            display(['AICNCC method has gotten TOF diff and used ',num2str(toc),' s']);

            
        % Cross correlation method for TOF difference
        elseif strcmp(varargin{5},'CC')
            display('start running cross correlation method for TOF difference');
            tic;
            tof_map = tof_picker_CC_diff(varargin{1}, varargin{2}, varargin{3}, varargin{4}, speed_in_water, ring_diameter);
            display(['cross correlation method has gotten TOF diff and used ',num2str(toc),' s']);
        else
            display(['there is no the method, ',varargin{5},', please checke the method name!!!']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First peak method (TH) for TOF map calculation
% 
% input: 
%   1. recieve_data: 3D received data. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   2. time_hist:  2D data. 1stD is time, 2nd is transmit elements.
%
% Output:
%   1. tof_map: TOF calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tof_map = tof_picker_FP(receive_data,time_hist)

%     [data_sample_num, receive_num, transmit_num] = size(receive_data);                    % get receve data size.
%     
%     [X, Y, Z] = ndgrid(1:data_sample_num, 1:receive_num, 1:transmit_num);
%     F = griddedInterpolant(X,Y,Z,receive_data);
%     [Xq, Yq, Zq] = ndgrid(1:0.5:data_sample_num, 1:receive_num, 1:transmit_num);
%     receive_data = F(Xq,Yq,Zq);
%     
%     
%     [X, Y] = ndgrid(1:data_sample_num, 1:transmit_num);
%     F = griddedInterpolant(X,Y,time_hist);
%     [Xq, Yq] = ndgrid(1:0.5:data_sample_num, 1:transmit_num);
%     time_hist = F(Xq,Yq);
%     clear Xq Yq Zq;


    % change time hist of water and pattern same. and do interpolation for higher accuracy.
    [data_sample_num, receive_num, transmit_num] = size(receive_data); 
    dt = time_hist(2) - time_hist(1); 
    time_hist_new = repmat((dt/2:dt/2:max(time_hist(:,1)))',[1,transmit_num]);
    receive_data = interp1(time_hist(:,1),receive_data, time_hist_new(:,1),'linear','extrap');               % receive water data interpolation
    
    time_hist = time_hist_new; clear time_hist_new dt;
    [data_sample_num, receive_num, transmit_num] = size(receive_data); 
    
    window_size = 16;
    
    % threshold method to get the start position for first peak method.
    receive_data = abs(receive_data);
    threshold_map = max(receive_data)./90;                                                                        % get threshold mat for each pair of receive and transmit.
    seg_map = receive_data > repmat(reshape(threshold_map,[1,receive_num,transmit_num]),[data_sample_num,1,1]);     % segmetate all receive data using threshld map.
    [m,first_arrival_map] = max(seg_map,[],1);                                                                                % get the position of pulse
    
    % get window data for first peak detection.
    window_data_offset_start = reshape(first_arrival_map,[receive_num,transmit_num]) + repmat((1:receive_num)'-1, [1, transmit_num]).*data_sample_num + repmat((1:transmit_num)-1,[receive_num,1]).*receive_num.*data_sample_num;
    window_data_offset = repmat(reshape(window_data_offset_start,[1,receive_num,transmit_num]), [window_size,1,1]) + repmat((1:window_size)'-1, [1, receive_num,transmit_num]);
    window_data = receive_data(window_data_offset);
    
    % allocate space for tof
    first_peak_arrival_map = zeros(receive_num,transmit_num);
    
    % fist peak detection.
    for transmit_i = 1:transmit_num
        for receive_i = 1:receive_num
%             [peak_value,peak_index] = findpeaks(window_data(:,receive_i,transmit_i)); %%% find first peak;
%             first_peak_arrival_map(receive_i,transmit_i) = peak_index(1);
            [peak_value,peak_index] = max(window_data(:,receive_i,transmit_i));   %%% find maximum
            first_peak_arrival_map(receive_i,transmit_i) = peak_index(1);
        end
    end
    
    % get tof.
    first_peak_arrival_map = first_peak_arrival_map +  reshape(first_arrival_map,[receive_num,transmit_num]) -1;
    tof_map =first_peak_arrival_map.* (time_hist(2,1) - time_hist(1,1));                      % calculate the tof using positon of pulse.
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Threshold method (TH) for TOF map calculation
% 
% input: 
%   1. recieve_data: 3D received data. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   2. time_hist:  2D data. 1stD is time, 2nd is transmit elements.
%
% Output:
%   1. tof_map: TOF calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tof_map = tof_picker_TH(receive_data, time_hist)

    [data_sample_num, receive_num, transmit_num] = size(receive_data);                                              % get receve data size.
    
    receive_data = abs(receive_data);
    threshold_map = max(receive_data)./90;                                                                        % get threshold mat for each pair of receive and transmit.
    seg_map = receive_data > repmat(reshape(threshold_map,[1,receive_num,transmit_num]),[data_sample_num,1,1]);     % segmetate all receive data using threshld map.
    [m,tof_map] = max(seg_map,[],1);                                                                                % get the position of pulse
    tof_map = reshape(tof_map,[receive_num,transmit_num]).* (time_hist(2,1) - time_hist(1,1));                      % calculate the tof using positon of pulse.
 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Akaike information criterion(AIC) method for TOF difference calculation 
% 
% input: 
%   1. recieve_data_pattern:    3D received data for pattern. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   2. time_hist_pattern:       2D data for pattern. 1stD is time, 2nd is transmit elements.
%   4. speed_in_water:          sound speed of water (m/s), during obtain receive_data_water and time_hist_water.
%   5. ring_diameter:           ring diameter (m) for ultrasound CT.
%
% Output:
%   1. tof_map: TOF calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tof_map = tof_picker_AIC(receive_data, time_hist, speed_in_water, ring_diameter)
       
    [data_sample_num, receive_num, transmit_num] = size(receive_data);      % get size of received data. 

    % claculate time of flight for water
    delta_angle = 2*pi/receive_num;                         % delta angle between neighbor elements
    rotation_angles = (0:receive_num-1).*delta_angle;       % rotation angle list for each elements
    element_pos_x = ring_diameter/2*cos(rotation_angles);                           % elements pos x, first element position was setted to be (ring_diameter/2, 0)
    element_pos_y = ring_diameter/2*sin(rotation_angles);                           % elements pos y 
    element_dis = sqrt((ring_diameter/2-element_pos_x).^2 + (0-element_pos_y).^2);  % element distance to the first elements
    tof_water = element_dis./speed_in_water;                                        % time of flight from first element to each elements. 

    % prepare parameters for tof calculation
    tof_map = zeros(receive_num,transmit_num);                      % allocate space for result, tof_map
    window_size_half = 128;                                         % AIC half window size. two times of it are same to search range for 10MHz
    AIC = ones(2*window_size_half,receive_num,transmit_num)*inf;    % AIC matrix initialization which is for saving calculated AIC value.
    
    
 %   AIC_temp = AIC;
    
    
    % loop for transmit elements
    for transmit_i = 1:transmit_num
        
        % read ch numbers line receive data for current transmition
        dt = time_hist(2,transmit_i) - time_hist(1,transmit_i);             % delta time 
        tof_water_rotated = circshift(tof_water', transmit_i-1)';           % circle shift for find right time for eache receive element of current transmition.
        [m_value, m_r] = max(tof_water_rotated);                            % find the position of latest arrival receive signal.
        
        curr_receive_data = receive_data(:,:,transmit_i);                   % get all receive data for current transmit element
        curr_receive_data_extent = cat(1, repmat(curr_receive_data(1:window_size_half,m_r),[1,receive_num]) , curr_receive_data, curr_receive_data(1:window_size_half,:));  % extent the receive data to avoid out of array bounds
        
        % calculated AIC window center and window range using tof of water. 
        window_center = round(tof_water_rotated/dt) + window_size_half;         % window centers for each receive element of current transmition (Position in extented receive data matrix).       
        window_start = window_center - window_size_half + 1;                    % window start position for each receive element
        
        % setting search region
        scan_region_start = round((tof_water_rotated/dt).*0.9) + window_size_half - window_start;   % search region start positon in window data
        scan_region_end = round((tof_water_rotated/dt).*1.1) + window_size_half - window_start+1;   % search region end position in window data 
        
        % obtain AIC window data for each receive element.
        window_data_offset_start = ((1:receive_num)-1)*size(curr_receive_data_extent,1) + window_start;                                     % window start postion in the extented receive data for each receive element   
        window_data_offset = repmat(window_data_offset_start,[2*window_size_half,1]) + repmat((0:2*window_size_half-1)', [1,receive_num]);  % all window data position in the extented receive data.
        window_data = curr_receive_data_extent(window_data_offset);                                                                         % get all window data
    
        % prepare parameters for AIC value calculation
        scan_region_mask = zeros(1,receive_num);    % initialize search region mask.
        var_lim_mask = zeros(1,receive_num);        % limite the var value to make sure var of 2nd segment larger than the 1st segment.
        
        % AIC value calculation
        for AIC_i = 1:2*window_size_half 
            
            scan_region_mask(:) = 0;                % initialize for each loop.
            scan_region_mask(find((AIC_i<scan_region_start) +(AIC_i> scan_region_end))) = inf;      % value of positon outof search region will be setted to be inf.
            
            if AIC_i == 2*window_size_half
                continue;
            end
            
            S1 = var(window_data(1:AIC_i,:));       % variance of head half window
            S2 = var(window_data(AIC_i+1:end,:));   % variance of tail half winwon
            
            
            S1(S1<1) = 1; % avoid negtive, since log 0 is -inf
            S2(S2<1) = 1;

            % limitation mask for var value comparison.
            var_lim_mask(:) = 0;
            var_lim_mask(find(S1>S2)) = inf;
            
            AIC(AIC_i,:,transmit_i) = AIC_i*log(S1) + (2*window_size_half-AIC_i-1)*log(S2) + scan_region_mask + var_lim_mask; % AIC value map and value of psotion outof search region will be setted to be inf.
 %           AIC_temp(AIC_i,:,transmit_i) = AIC_i*log(S1) + (2*window_size_half-AIC_i-1)*log(S2); 
        end

        [AIC_min,p] = min(AIC(:,:,transmit_i));                                     % find min AIC value position in window for each receive element
        tof_map(:,transmit_i) = ((window_start + p - 1) - window_size_half)*dt ;    %(window_start + p - 1)is result in extented coordinate, to change it to original coordinate,  - window_size_half.
    
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted akaike information criterion (WAIC) method for TOF difference calculation 
%   Note: this method was proposed by Duric group, pubulished in "Ultrasonics Vol. 49, Iss. 1, Jan. 2009, P. 61-72" 
%
% input: 
%   1. recieve_data_pattern:    3D received data for pattern. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   2. time_hist_pattern:       2D data for pattern. 1stD is time, 2nd is transmit elements.
%   4. speed_in_water:          sound speed of water (m/s), during obtain receive_data_water and time_hist_water.
%   5. ring_diameter:           ring diameter (m) for ultrasound CT.
%
% Output:
%   1. tof_map: TOF calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tof_map = tof_picker_WAIC(receive_data, time_hist, speed_in_water, ring_diameter)

    [data_sample_num, receive_num, transmit_num] = size(receive_data);      % get size of received data. 

    % claculate time of flight for water
    delta_angle = 2*pi/receive_num;                         % delta angle between neighbor elements
    rotation_angles = (0:receive_num-1).*delta_angle;       % rotation angle list for each elements
    element_pos_x = ring_diameter/2*cos(rotation_angles);                           % elements pos x, first element position was setted to be (ring_diameter/2, 0)
    element_pos_y = ring_diameter/2*sin(rotation_angles);                           % elements pos y 
    element_dis = sqrt((ring_diameter/2-element_pos_x).^2 + (0-element_pos_y).^2);  % element distance to the first elements
    tof_water = element_dis./speed_in_water;                                        % time of flight from first element to each elements. 

    % prepare parameters for tof calculation
    tof_map = zeros(receive_num,transmit_num);                      % allocate space for result, tof_map
    window_size_half = 100;                                         % AIC half window size. two times of it are same to search range.
    AIC = ones(2*window_size_half,receive_num,transmit_num)*inf;   % AIC matrix initialization which is for saving calculated AIC value.

    % loop for transmit elements
    for transmit_i = 1:transmit_num
        
        % read ch numbers line receive data for current transmition
        dt = time_hist(2,transmit_i) - time_hist(1,transmit_i);             % delta time 
        tof_water_rotated = circshift(tof_water', transmit_i-1)';           % circle shift for find right time for eache receive element of current transmition.
        [m_value, m_r] = max(tof_water_rotated);                              % find the position of latest arrival receive signal.
        
        curr_receive_data = receive_data(:,:,transmit_i);                           % get all receive data for current transmit element
        curr_receive_data_extent = cat(1, repmat(curr_receive_data(1:window_size_half,m_r),[1,receive_num]) , curr_receive_data, curr_receive_data(1:window_size_half,:));  % extent the receive data to avoid out of array bounds
        
        % calculated AIC window center and window range using tof of water. 
        window_center = round(tof_water_rotated/dt) + window_size_half;       % window centers for each receive element of current transmition (Position in extented receive data matrix).       
        window_start = window_center - window_size_half + 1;                % window start position for each receive element
        
        % setting search region
        scan_region_start = round((tof_water_rotated/dt).*0.9) + window_size_half - window_start;   % search region start positon in window data
        scan_region_end = round((tof_water_rotated/dt).*1.1) + window_size_half - window_start+1;   % search region end position in window data 
        
    
        % obtain AIC window data for each receive element.
        window_data_offset_start = ((1:receive_num)-1)*size(curr_receive_data_extent,1) + window_start;                                     % window start postion in the extented receive data for each receive element   
        window_data_offset = repmat(window_data_offset_start,[2*window_size_half,1]) + repmat((0:2*window_size_half-1)', [1,receive_num]);  % all window data position in the extented receive data.
        window_data = curr_receive_data_extent(window_data_offset);                                                                         % get all window data
    
        % prepare parameters for AIC value calculation
        scan_region_mask = zeros(1,receive_num);    % initialize search region mask.
        
        % AIC value calculation
        for AIC_i = 1:2*window_size_half 
            
            scan_region_mask(:) = 0;                % initialize for each loop.
            scan_region_mask(find((AIC_i<scan_region_start) +(AIC_i> scan_region_end))) = inf;      % value of positon outof search region will be setted to be inf.
            
            if AIC_i == 2*window_size_half
                continue;
            end
            
            S1 = var(window_data(1:AIC_i,:));       % variance of head half window
            S2 = var(window_data(AIC_i+1:end,:));   % variance of tail half winwon
        
            S1(S1<1) = 1; % avoid negtive, since log 0 is -inf
            S2(S2<1) = 1;
        
            AIC(AIC_i,:,transmit_i) = AIC_i*log(S1) + (2*window_size_half-AIC_i-1)*log(S2) + scan_region_mask; % AIC value map and value of psotion outof search region will be setted to be inf.
        end
    
        AIC_min = min(AIC(:,:,transmit_i));
        
        % calculate the weight for each AIC value
        exp_delta_AIC = exp(-(AIC(:,:,transmit_i)- repmat(AIC_min,[2*window_size_half,1])));    % find min AIC value in window for each receive element
        weight_normal = exp_delta_AIC./repmat(sum(exp_delta_AIC,1),[2*window_size_half,1]);     % calculate weight for each AIC value;
        p = sum(weight_normal .* repmat((1:2*window_size_half)',[1 receive_num]));              % calculated weighted min AIC position.
        tof_map(:,transmit_i) = ((window_start + p - 1) - window_size_half)*dt ;                %(window_start + p - 1)is result in extented coordinate, to change it to original coordinate, - window_size_half.
    
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algorithm combined by akaike information criterion and neighbor cross correlation method and named by (AICCC) for TOF difference calculation 
%   Note: this method was proposed by our team in June 2014.
%
% input: 
%   1. recieve_data_pattern:    3D received data for pattern. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   2. time_hist_pattern:       2D data for pattern. 1stD is time, 2nd is transmit elements.
%   4. speed_in_water:          sound speed of water (m/s), during obtain receive_data_water and time_hist_water.
%   5. ring_diameter:           ring diameter (m) for ultrasound CT.
%
% Output:
%   1. tof_map: TOF calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tof_map = tof_picker_AICNCC(receive_data, time_hist, speed_in_water, ring_diameter)
   
    % prepare input parematers for AICC
    [data_sample_num, receive_num, transmit_num] = size(receive_data);                          % get size of received data. 
    tof_map_AIC = tof_picker_AIC(receive_data, time_hist, speed_in_water, ring_diameter);       % get tof map calculated by AIC method.
    dt = time_hist(2,1) - time_hist(1,1);                                                       % get delta t
    window_size = 64;                                                                           % cross correlation window size
    scan_region_size =30;                                                                       % cross correlation search region size.    % Note: in my understand, this value should be half T lenght of receive puls 19th Jun. 2014.
    receive_data_extent = cat(1, receive_data, zeros(window_size, receive_num, transmit_num));  % extent receive data lenght for avoiding array outof boundary
    
    % prepare cross correlation window data.
    window_start = round(tof_map_AIC./dt);          % windnow start position. which is same to the AIC results positon
    window_start(window_start<1) = 1;               % avoiding array outof boundary.
    window_start_offset = window_start + repmat((0:receive_num-1)',[1,transmit_num]).* (data_sample_num+window_size) + repmat(0:transmit_num-1,[receive_num,1]) .* (data_sample_num+window_size) .* receive_num;    % window start offset in extented receive data for each pair of tranmit and receive.
    window_data = receive_data_extent(repmat(reshape(window_start_offset, 1, receive_num, transmit_num),[window_size,1,1]) + repmat((0:window_size-1)',[1, receive_num, transmit_num]));                            % obtaining window data from extented receive data.
    
    % prepare cross correlation search region data
    scan_region_start = window_start - round(scan_region_size./2) + 1;      % search region start position
    scan_region_end = window_start +  round(scan_region_size./2);           % search region end position
    scan_region_start (scan_region_start<1) = 1;                            % avoiding array outof boundary.
    scan_region_start_offset = scan_region_start + repmat((0:receive_num-1)',[1,transmit_num]).* (data_sample_num+window_size) + repmat(0:transmit_num-1,[receive_num,1]) .* (data_sample_num+window_size) .* receive_num;                      % search region start offset in extented receive data for each pair of tranmit and receive.
    scan_region_data = receive_data_extent(repmat(reshape(scan_region_start_offset, 1, receive_num, transmit_num),[scan_region_size+2+window_size,1,1])   +   repmat((0:scan_region_size+2+window_size-1)',[1, receive_num, transmit_num]));    % obtaining search region data from extented receive data.
    
    % cross correlation neighbor region parameter setting
    neighbor_range = 3;                                 % neighbor region side lenght. Note neighbor region is a square rectangle.
    neighbor_num = neighbor_range.^2;                   % neighbor region size.
    neighbor_center = round((neighbor_range+1)/2);      % neighbor region center position. which is in the center of the square rectangle and is locted object data.
    p_shift = zeros(neighbor_num, receive_num, transmit_num);                           % allocate space for shift results inspected by each neighbor cross correlation. 
    tof_neighbor_CC = zeros(neighbor_num, receive_num, transmit_num);                   % allocate space for tof inspected by each neighbor cross correlation.
    cross_correlation_max_map = zeros(neighbor_num, receive_num, transmit_num);         % allocate space for maximum correlation inspected by each neighbor cross correlation.
    curr_cross_correlation_map = zeros(scan_region_size, receive_num, transmit_num);    % allocate spaace for cross correlation results for same neighbor's window data in differen position of scan region data.
    
    % loop for finding each neighbor's shift and maximum correlation by cross correlation.
    for neighbor_receive_i = 1:neighbor_range
        for neighbor_transmit_i = 1:neighbor_range
            
            curr_neighbor_counter = neighbor_receive_i + (neighbor_transmit_i-1).*neighbor_range;   % current neighbor counter number (number code in 1D, for example, 1-9 for 3*3 neighbor region)
            
            % if current neighbor region is center position, which is not neighbor but itself.
            if neighbor_receive_i == neighbor_center && neighbor_transmit_i == neighbor_center
                cross_correlation_max_map(curr_neighbor_counter,:,:) = 1;       % cross correlation should be 1, since it is itself
                tof_neighbor_CC(curr_neighbor_counter,:,:) = tof_map_AIC;               % tof shift should be AIC tof map
                continue;
            end
            
            % prepare data for cross correlation map calculation of current neighbor
            window_data_shift = circshift(window_data, [0, neighbor_receive_i-neighbor_center, neighbor_transmit_i-neighbor_center]);   % shift current neighbor window data to make it corresponding to the center search region. 
            curr_cross_correlation_map(:) = 0;      % initialize current_cross_correlation map with 0;
            
            % loop for calculate the cross correlation map, which is correlation result of current neighbor window in different position of center search region..
            for ii = 1 : scan_region_size
                scan_region_mask = scan_region_end > ii;    % end mask. check the search region is not ended
                
                % normalized cross correlation calculation.
                curr_cross_correlation_map(ii,:,:) = sum(window_data_shift.*scan_region_data(ii:ii+window_size-1,:,:))./sqrt(sum(window_data_shift.^2).*sum(scan_region_data(ii:ii+window_size-1,:,:).^2)).* reshape(scan_region_mask,[1,receive_num,transmit_num]);     
            end
            
            curr_cross_correlation_map (curr_cross_correlation_map<=0) = 0; % chang too low correlation to be zero.   this is an important threshold to remove affection of cross correlation match failed neighbor.
            curr_cross_correlation_map = curr_cross_correlation_map.^10;    % testing
            
            [cross_correlation_max_map(curr_neighbor_counter,:,:), p_shift(curr_neighbor_counter,:,:)] = max (curr_cross_correlation_map);                      % get maximum correlation and shift position of current neighbor.
            tof_neighbor_CC(curr_neighbor_counter,:,:) = (reshape(p_shift(curr_neighbor_counter,:,:),[receive_num, transmit_num]) + scan_region_start-1).*dt;   % get time of filght calculated by crrent neighbor cross correlation.
            
        end
    end
    
    % get last time of flight result by each tof calucated by each neighbor with similarity weighted.
    tof_map = reshape(sum(tof_neighbor_CC.*(cross_correlation_max_map./repmat(sum(cross_correlation_max_map),[neighbor_num,1,1]))),[receive_num, transmit_num]);
   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross Correlation method for TOF difference calculation
%   Note: The cross correlation between received water data and received pattern data will be calculated. 
%   
%   History:
%   1. original wittern by Xiaolei Qu in 2014 to fast implement nakamura method.
%   2. modified by Xiaolei Qu in 2015 to move window center to water arrival time.
% 
% input: 
%   1. recieve_data_pattern:    3D received data for pattern. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   2. time_hist_pattern:       2D data for pattern. 1stD is time, 2nd is transmit elements.
%   3. receive_data_water:      3D received data for water. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
%   4. time_hist_water:         2D data for water. 1stD is time, 2nd is transmit elements.
%   5. speed_in_water:          sound speed of water (m/s), during obtain receive_data_water and time_hist_water.
%   6. ring_diameter:           ring diameter (m) for ultrasound CT.
%
% Output:
%   1. tof_dif_map: TOF difference calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tof_dif_map = tof_picker_CC_diff(receive_data_pattern, time_hist_pattern, receive_data_water, time_hist_water, speed_in_water, ring_diameter)

%     % change time hist of water and pattern same. and do interpolation for higher accuracy.
%     [data_sample_num, receive_num, transmit_num] = size(receive_data_pattern);                      % get size of receve data matrix
%     time_hist = linspace(min(min(time_hist_pattern(:)),min(time_hist_water(:))),max(max(time_hist_pattern(:)),max(time_hist_water(:))),2*data_sample_num);  % get the interpolated time hist
%     %time_hist= linspace(time_hist_water(1,transmit_num),time_hist_water(end,transmit_num),2*size(receive_data_water,1));
%     receive_data_water = interp1(time_hist_water(:,1),receive_data_water, time_hist,'linear','extrap');               % receive water data interpolation
%     receive_data_pattern = interp1(time_hist_pattern(:,1),receive_data_pattern, time_hist,'linear','extrap');         % receive pattern data interpolation
%     dt = time_hist(2) - time_hist(1);                                                               % delta t

    time_hist = time_hist_water(:,1);
%     receive_data_pattern = interp1(time_hist_pattern(:,1),receive_data_pattern, time_hist,'linear','extrap');
    [data_sample_num, receive_num, transmit_num] = size(receive_data_pattern);
    dt = time_hist(2) - time_hist(1);
    
    % calculated idea arrival time for water
    element_num = receive_num;                              % receive element numbre of the Ring
    delta_angle = 2*pi/element_num;                         % delta angle between neighbor elements
    rotation_angles = (0:element_num-1).*delta_angle;       % rotation angle list for each elements
    element_pos_x = ring_diameter/2*cos(rotation_angles);                           % elements pos x
    element_pos_y = ring_diameter/2*sin(rotation_angles);                           % elements pos y 
    element_dis = sqrt((ring_diameter/2-element_pos_x).^2 + (0-element_pos_y).^2);  % element distance to the first elements
    tof_water = element_dis./speed_in_water;                                        % time of flight from first element to each elements. 
    
    % parameter prepare for first break picking
    tof_dif_map = zeros(receive_num,transmit_num);              % allocate space for time of flight map which will be last result of this function.
    max_similarity_map = zeros(receive_num, transmit_num);      % allocate space for maximum similarity map.
    window_size = 128;                                          % cross correlation window size.
    scan_region_size = round((max(tof_water) *0.2)/dt);         % scan region size or search region size. pattern sond speed range is water speed + or - 10%.
    correlation_map = zeros(scan_region_size, receive_num);     % allocate space for correlation map 
    
    % loop for find time of flight map
    for transmit_i = 1:transmit_num
        
        curr_receive_data_pattern = receive_data_pattern(:,:,transmit_i);   % get pattern data for current transmiter element
        curr_receive_data_water = receive_data_water(:,:,transmit_i);       % get water data for current transmiter element
        curr_receive_data_pattern_extent = cat(1, zeros(window_size,receive_num), curr_receive_data_pattern, zeros(window_size,receive_num));   % extent pattern data to avoid array outof bounds 
        curr_receive_data_water_extent = cat(1, zeros(window_size,receive_num), curr_receive_data_water, zeros(window_size,receive_num));       % extent water data to avoid array outof bounds
    
        
        % window data from received water data for cross correlation 
        tof_water_rotated = circshift(tof_water', transmit_i-1)'+ window_size.*dt;      % arrival time in water for each receive element
        window_start_water = round(tof_water_rotated/dt)-round(window_size/2);          % window start position was changed to make sure window center in water arrival time. 
        window_start_water(window_start_water<1) = 1;                                   %for avoiding zero index which only happen when the transmit and receive elements are same one.
        window_start_water_offset = ((1:receive_num)-1)*size(curr_receive_data_water_extent,1) + window_start_water;    % start position of each recieve elments in 1D expression.
        window_data_water = curr_receive_data_water_extent(repmat(window_start_water_offset,[window_size,1]) + repmat((0:window_size-1)', [1,receive_num]));    % window data obtained for each receive.
        
    
        % scanning region data from receive pattern data.
        scan_region_start_pattern = round(tof_water_rotated/dt.*0.9)-round(window_size/2);                             % start position. pattern sound speed must larger than 90% water sound speed.
        scan_region_end_pattern = round(tof_water_rotated/dt*1.1) - scan_region_start_pattern;    % end position
        scan_region_start_pattern (scan_region_start_pattern<1) = 1;                            % for avoiding zero index.
        scan_region_start_pattern_offset = ((1:receive_num)-1).*size(curr_receive_data_pattern_extent,1) + scan_region_start_pattern;   % start position in 1D expression
        scan_region_data_pattern = curr_receive_data_pattern_extent(repmat(scan_region_start_pattern_offset,[scan_region_size+window_size,1]) + repmat((0:(scan_region_size+window_size-1))', [1,receive_num]));    % scan region data obtaining
        
        
        % correlation_map initialization by 0;
        correlation_map(:) = 0;
        % loop for correlation mat calculation.
        for ii = 1:scan_region_size
            scan_region_mask = scan_region_end_pattern > ii;    % mask for removing receive elements outof search region.
            correlation_map(ii,:) = (sum((window_data_water.*scan_region_data_pattern(ii:(ii+window_size-1),:)),1)./sqrt(abs(sum(window_data_water.^2,1).*sum(scan_region_data_pattern(ii:(ii+window_size-1),:).^2,1)))).*scan_region_mask;
        end
        
         
        % get time of flight differece map between water and pattern
        [m, p] = max(correlation_map,[],1);     % get maximum correlation position for each receive elements
        tof_dif_map(:,transmit_i) = ((p + scan_region_start_pattern - 1) - window_start_water).*dt;     % get time of flight difference for each recieve elements.
        max_similarity_map(:,transmit_i) = m;   % maximimu similarity for each time different calculation.
    
    end

end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Cross Correlation method for TOF difference calculation
% %   Note: The cross correlation between received water data and received pattern data will be calculated. 
% % 
% % input: 
% %   1. recieve_data_pattern:    3D received data for pattern. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
% %   2. time_hist_pattern:       2D data for pattern. 1stD is time, 2nd is transmit elements.
% %   3. receive_data_water:      3D received data for water. 1stD is time, 2nd is receive elements, 3rd is transmit elements.
% %   4. time_hist_water:         2D data for water. 1stD is time, 2nd is transmit elements.
% %   5. speed_in_water:          sound speed of water (m/s), during obtain receive_data_water and time_hist_water.
% %   6. ring_diameter:           ring diameter (m) for ultrasound CT.
% %
% % Output:
% %   1. tof_dif_map: TOF difference calculation result. 2D data. 1st D is recieve elements, 2nd is transmit elements.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function tof_dif_map = tof_picker_CC_diff(receive_data_pattern, time_hist_pattern, receive_data_water, time_hist_water, speed_in_water, ring_diameter)
% 
%     % change time hist of water and pattern same. and do interpolation for higher accuracy.
%     [data_sample_num, receive_num, transmit_num] = size(receive_data_pattern);                      % get size of receve data matrix
%     time_hist = linspace(min(min(time_hist_pattern(:)),min(time_hist_water(:))),max(max(time_hist_pattern(:)),max(time_hist_water(:))),2*data_sample_num);  % get the interpolated time hist
%     %time_hist= linspace(time_hist_water(1,transmit_num),time_hist_water(end,transmit_num),2*size(receive_data_water,1));
%     receive_data_water = interp1(time_hist_water(:,1),receive_data_water, time_hist);               % receive water data interpolation
%     receive_data_pattern = interp1(time_hist_pattern(:,1),receive_data_pattern, time_hist);         % receive pattern data interpolation
%     dt = time_hist(2) - time_hist(1);                                                               % delta t
%     
%     
%     % calculated idea arrival time for water
% %     speed_in_water = 1523;                              % m/s sound speed of water, for reference window center.
% %     ring_diameter = 0.035;                              % 3.5cm diameter of Ring
%     element_num = receive_num;                          % receive element numbre of the Ring
%     delta_angle = 2*pi/element_num;                     % delta angle between neighbor elements
%     rotation_angles = (0:element_num-1).*delta_angle;   % rotation angle list for each elements
%     element_pos_x = ring_diameter/2*cos(rotation_angles);                           % elements pos x
%     element_pos_y = ring_diameter/2*sin(rotation_angles);                           % elements pos y 
%     element_dis = sqrt((ring_diameter/2-element_pos_x).^2 + (0-element_pos_y).^2);  % element distance to the first elements
%     tof_water = element_dis./speed_in_water;                                        % time of flight from first element to each elements. 
%     
%     % parameter prepare for first break picking
%     tof_dif_map = zeros(receive_num,transmit_num);              % allocate space for time of flight map which will be last result of this function.
%     window_size = 128;                                          % cross correlation window size.
%     scan_region_size = round((max(tof_water) *0.2)/dt);         % scan region size or search region size. pattern sond speed range is water speed + or - 10%.
%     correlation_map = zeros(scan_region_size, receive_num);     % allocate space for correlation map 
%     
%     % loop for find time of flight map
%     for transmit_i = 1:transmit_num
%         
%         curr_receive_data_pattern = receive_data_pattern(:,:,transmit_i);   % get pattern data for current transmiter element
%         curr_receive_data_water = receive_data_water(:,:,transmit_i);       % get water data for current transmiter element
%         curr_receive_data_pattern_extent = cat(1, curr_receive_data_pattern, zeros(window_size,receive_num));   % extent pattern data to avoid array outof bounds 
%         curr_receive_data_water_extent = cat(1, curr_receive_data_water, zeros(window_size,receive_num));       % extent water data to avoid array outof bounds
%     
%         
%         % window data from received water data for cross correlation 
%         tof_water_rotated = circshift(tof_water', transmit_i-1)';           % arrival time in water for each receive element
%         window_start_water = round(tof_water_rotated/dt);                   % start position for arrival in water.
%         window_start_water(window_start_water<1) = 1;                       %for avoiding zero index which only happen when the transmit and receive elements are same one.
%         window_start_water_offset = ((1:receive_num)-1)*size(curr_receive_data_water_extent,1) + window_start_water;    % start position of each recieve elments in 1D expression.
%         window_data_water = curr_receive_data_water_extent(repmat(window_start_water_offset,[window_size,1]) + repmat((0:window_size-1)', [1,receive_num]));    % window data obtained for each receive.
%         
%     
%         % scanning region data from receive pattern data.
%         scan_region_start_pattern = round(window_start_water.*0.9);                             % start position. pattern sound speed must larger than 90% water sound speed.
%         scan_region_end_pattern = round(window_start_water*1.1) - scan_region_start_pattern;    % end position
%         scan_region_start_pattern (scan_region_start_pattern<1) = 1;                            % for avoiding zero index.
%         scan_region_start_pattern_offset = ((1:receive_num)-1).*size(curr_receive_data_pattern_extent,1) + scan_region_start_pattern;   % start position in 1D expression
%         scan_region_data_pattern = curr_receive_data_pattern_extent(repmat(scan_region_start_pattern_offset,[scan_region_size+window_size,1]) + repmat((0:(scan_region_size+window_size-1))', [1,receive_num]));    % scan region data obtaining
%         
%         
%         % correlation_map initialization by 0;
%         correlation_map(:) = 0;
%         % loop for correlation mat calculation.
%         for ii = 1:scan_region_size
%             scan_region_mask = scan_region_end_pattern > ii;    % mask for removing receive elements outof search region.
%             correlation_map(ii,:) = sum((window_data_water.*scan_region_data_pattern(ii:(ii+window_size-1),:)),1).*scan_region_mask;
%         end
%         
%          
%         % get time of flight differece map between water and pattern
%         [m, p] = max(correlation_map,[],1);     % get maximum correlation position for each receive elements
%         tof_dif_map(:,transmit_i) = ((p + scan_region_start_pattern - 1) - window_start_water).*dt;     % get time of flight difference for each recieve elements.
%         
%     end
% 
% end









