
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          image reconstruction
%   Description:    ultrasound sound speed tomography image reconstruction straight ray method and ray tracing method. 
%   Author:         Xiaolei Qu
%   Date:           2014.08.06
%   Version:        1.0
%   History:        1.0  straight ray (S) method and bent ray brute force method (BBF) is writtern basing on code of nakamura san.
%
%   Input:
%       pjt_data        travel time differece map. 2D data 
%                       1st D. is receiver, 
%                       2nd D. is emitter, 
%       method_type     method_type:    
%                       S: straight ray reconstruction
%                       CBF: VERY SLOW METHOD: curve(bent) ray reconstruction method, using brute force method for ray linking.
%                       CNR: curve(bent) ray reconstruction method, using Newton-Raphson method for ray linking.
%                       CBS: curve(bent) ray reconstruction method, using Bisection method for ray linking.
%                       CVR: curve(bent) ray reconstruction method, using virtual receive technique avoding ray linking.
% 
%   Output:
%       recon_image:    reconstructed image (2D slice)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recon_image= img_reconstruct(varargin)


    if nargin == 2
        pjt_data = varargin{1};
        method_type = varargin{2};
        
        iter_result_path = pwd;             % iteration result saving path.
        
    elseif nargin == 3
        pjt_data = varargin{1};
        method_type = varargin{2};
        iter_result_path = varargin{3};
        
    else
        display('input parameter not enough');
        recon_image = -1;
        return;
    end
        
    
    

    % prepare important parameters for recon
    speed_in_water = 1520;          % 1540 for simulation data, 1482 for 20141014experiment data.
    ring_diameter = 0.2;            % ring diameter
    image_size_m = 0.21;            % size of image in meters
    
    element_num = 256;              % element numbers
    image_size_p = 256;             % recon image size in pixel. 20160630ˆäãC³@256¨1024
    
    % angle limited in reconstruction
    available_receiver_angle_range = 5/3*pi;
    
    % calculate iteration times for straight ray reconstruction
    iteration_num_s_calculation = 20;
     
    % iteration reconstruction parameters for curve ray reconstruction
    iteration_num_c_s_calculation = 5;          % straight ray calculation iteration times for depature points of curve ray tracing.
    iteration_num_c_calculation = 2;            % curve ray calculation iteration times for each ray linking results.
    iteration_num_c_ray_tracing = 100;           % curve ray linking iteration times.
    
    % ray grid step size, generally it should be half of image resolution.
    ds = 2*speed_in_water/1.e6/4;       % too small steps will cause ray tracing problem. 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for removeing failed pjt data.
    % max difference of pjt for reciprocal emitter-receiver pair
    max_diff_pjt_emit_recei_pair = 3e-7;
    
    % the minimum reconstruction speed and maximum reconstruction speed for removing failed pjt data.
    min_recon_speed = 1200;
    max_recon_speed = 1650; %20160630ˆäãC³ 1600¨1650
    
    
    

    % straight ray reconstruction 
    if strcmp(method_type,'S')  
        
        display('start to calculate coefficients projection matrix ...');
        tic;
        [a, w] = straight_ray_tracing(ring_diameter, element_num, image_size_p, image_size_m);      % get coefficient projection matrix by straight ray tracing.
        display(['projection matrix calculation used ',num2str(toc),'s.']);
        display('start SART recon ...');
        recon_image = get_recon_image_SART(pjt_data, a, w, ones(image_size_p)*speed_in_water, image_size_p, image_size_m,iteration_num_s_calculation,1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed);    % reconstruction iteration calculation
        display(['S-SART recon image was obtained in ',num2str(toc),'s']);
    
    % curve ray recon    
    elseif strcmp(method_type,'CBF')  || strcmp(method_type,'CNR') || strcmp(method_type,'CBS') || strcmp(method_type,'CVR') || strcmp(method_type,'CFMM')
        
        % straight ray reconstruction for a departure point of curve ray reconstuction.
        display('start S-SART recon for initial depature point ...');
        tic;
        [a, w] = straight_ray_tracing(ring_diameter, element_num, image_size_p, image_size_m);                                                  % get coefficient projection matric by straight ray tracing
        %ˆäãC³ones(element_num)*speed_in_water¨ones(image_size_p)*speed_in_water
        recon_image = get_recon_image_SART(pjt_data, a, w, ones(image_size_p)*speed_in_water, image_size_p, image_size_m,iteration_num_c_s_calculation,1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed); % get recon image by S-SART method
        display(['S-SART recon image obtained in ',num2str(toc),'s']);
        
        % removing reconstruction noise closing to the transudcer ring.
        tempx= repmat((1:image_size_p)',[1,image_size_p]);
        tempy = repmat(1:image_size_p,[image_size_p,1]);
        mask_recon = sqrt((tempx-image_size_p/2).^2+(tempy-image_size_p/2).^2)>0.034/2/0.04*image_size_p;
        recon_image_median_filted = medfilt2(recon_image, [5 5]);
        recon_image(mask_recon) = recon_image_median_filted(mask_recon);
        
        display(['start C-SART-', method_type ,' recon ...']);
        for i=1:iteration_num_c_ray_tracing        % curve ray tracing iteration reconstruction,
            
            % Curve ray reconstruction using brute Force ray linking method
            if strcmp(method_type,'CBF')             
                [a, w] = curve_ray_tracing_CBF(ring_diameter, element_num, image_size_p, image_size_m, recon_image,ds,speed_in_water);                            % get coefficient projection matrix by curve ray tracing
                recon_image = get_recon_image_SART(pjt_data, a, w, recon_image,image_size_p, image_size_m,iteration_num_c_calculation, 1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed);      % get recon image by C-SART method
            
            % Curve ray reconstruction using Newton Raphson ray linking method
            elseif strcmp(method_type,'CNR')                
                [a, w] = curve_ray_tracing_CNR(ring_diameter, element_num, image_size_p, image_size_m, recon_image,ds,speed_in_water);                            % get coefficient projection matrix by curve ray tracing
                recon_image = get_recon_image_SART(pjt_data, a, w, recon_image,image_size_p, image_size_m,iteration_num_c_calculation, 1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed);      % get recon image by C-SART method
                
            % Curve ray reconstruction using bisection ray linking method
            elseif strcmp(method_type,'CBS')                
                [a, w] = curve_ray_tracing_CBS(ring_diameter, element_num, image_size_p, image_size_m, recon_image,ds,speed_in_water);                            % get coefficient projection matrix by curve ray tracing
                recon_image = get_recon_image_SART(pjt_data, a, w, recon_image,image_size_p, image_size_m,iteration_num_c_calculation, 1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed);      % get recon image by C-SART method
                
            % Curve ray reconstruction using Virtal Receiver ray linking method 
            elseif strcmp(method_type,'CVR')                
                [a, w, pjt_data_vr] = curve_ray_tracing_CVR(ring_diameter, element_num, image_size_p, image_size_m, recon_image ,pjt_data,ds,speed_in_water);     % get coefficient projection matrix by curve ray tracing
                recon_image = get_recon_image_SART(pjt_data_vr, a, w, recon_image,image_size_p, image_size_m,iteration_num_c_calculation, 1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed);  	% get recon image by C-SART method
                
            % wave front reconstruction using fast marching method (FMM). this method also is basing on ray theory.
            elseif strcmp(method_type,'CFMM')
                [a, w] = curve_wave_front_FMM(ring_diameter, element_num, image_size_p, image_size_m, recon_image,ds);                              % get coefficient projection matrix by curve ray tracing
                recon_image = get_recon_image_SART(pjt_data, a, w, recon_image,image_size_p, image_size_m,iteration_num_c_calculation, 1,available_receiver_angle_range,speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed);      % get recon image by C-SART method
                
            end
            
            % show current iteration time information.
            display([num2str(i),' times ray tracing overed in ',num2str(toc),'s' ]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing show and output current reconstruced images
%             figure; imagesc(recon_image,[1430 1650]); colorbar;
            save([iter_result_path,'\Method_',method_type,'_Iteration_',num2str(i),'.mat'],'recon_image');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
            
            % removing noise on the ring.
            recon_image_median_filted = medfilt2(recon_image, [5 5]);
            recon_image(mask_recon) = recon_image_median_filted(mask_recon);

        end
        display(['C-SART-',method_type,' recon image was obtained in',num2str(toc),'s.']);
        
    % method type input error    
    else
        % method type input mistake.
        recon_image = 0;
        display('method type input error');
    end 
    
end





function [a_ray_matrix, w_ray_matrix] = curve_wave_front_FMM(ring_diameter, element_num, image_size_p, image_size_m, curr_image,ds)

    image_resolution = image_size_m/image_size_p;
    ray_grid_step = ds;
    ray_grid_step = image_resolution/2;
    ray_grid_step = 0.5;
    
    
%     %%%%%%%%%%%%%%%% testing for smooth current image;
%     h = fspecial('average', [3 3]);     % average filter kernel
%     curr_image = filter2(h,curr_image);                   % average smoothed refractive index map.
    
    
    
    % calculate each element position, angle of start position is 0;
    theta_rg = linspace(0, 2*pi*(element_num-1)/element_num, element_num);
    ring_radius = ring_diameter/2;                      % radius of the ring.        
    element_pos_x = ring_radius * cos(theta_rg);        % element position x
    element_pos_y = ring_radius * sin(theta_rg);        % elmeent position y.

    element_pos_x_in_image = (image_size_p-1)/image_size_m * (element_pos_x+image_size_m./2) + 1;
    element_pos_y_in_image = (image_size_p-1)/image_size_m * (element_pos_y+image_size_m./2) + 1;
   
    
    ray_grid_pos_x = ones(1000,element_num,element_num).*inf;
    ray_grid_pos_y = ones(1000,element_num,element_num).*inf;
    
    
    for transmit_i = 1:element_num
        
        start_points = [round(element_pos_x_in_image(transmit_i)); round(element_pos_y_in_image(transmit_i))];
        
        D = msfm2d(curr_image./image_resolution,start_points,true,true);
        path = compute_path(D, [element_pos_x_in_image; element_pos_y_in_image ], ray_grid_step);
        
        
        for receive_i = 1:element_num
            N = size(path{receive_i},2);
            if N>1000
                N=1000;
                display('ray step too many');
            end
            
            % coordinate system transform.
            ray_grid_pos_x(1:N,receive_i,transmit_i) = ((path{receive_i}(1,1:N))'-1) * (image_size_m/(image_size_p-1)) - image_size_m/2;
            ray_grid_pos_y(1:N,receive_i,transmit_i) = ((path{receive_i}(2,1:N))'-1) * (image_size_m/(image_size_p-1)) - image_size_m/2;
        end
    
    end
    
    
    % get Coefficient matrix by curve ray tracing result.
    [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m);


end


function [T,Y]=msfm2d(F, SourcePoints, usesecond, usecross)

    % Distance image, also used to store the index of narrowband pixels 
    % during marching process
    T = zeros(size(F))-1;

    % Augmented Fast Marching (For skeletonize)
    Ed=nargout>1;

    % Euclidian distance image 
    if(Ed), Y = zeros(size(F)); end

    % Pixels which are processed and have a final distance are frozen
    Frozen   = zeros(size(F));

    % Free memory to store neighbours of the (segmented) region
    neg_free = 100000;
    neg_pos=0;
    if(Ed),
        neg_list = zeros(4,neg_free);
    else
        neg_list = zeros(3,neg_free);
    end

    % (There are 3 pixel classes:
    %   - frozen (processed)
    %   - narrow band (boundary) (in list to check for the next pixel with smallest distance)
    %   - far (not yet used)

    % Neighbours
    ne =[-1 0;
        1 0;
        0 -1;
        0 1];

    SourcePoints=int32(floor(SourcePoints));

    % set all starting points to distance zero and frozen
    for z=1:size(SourcePoints,2)
        % starting point
        x= SourcePoints(1,z); y=SourcePoints(2,z);
        % Set starting point to frozen and distance to zero
        Frozen(x,y)=1; T(x,y)=0;
    end

    % Add all neighbours of the starting points to narrow list
    for z=1:size(SourcePoints,2)
        % starting point
        x=SourcePoints(1,z); 
        y=SourcePoints(2,z);
        for k=1:4,
            % Location of neighbour
            i=x+ne(k,1); j=y+ne(k,2);
            % Check if current neighbour is not yet frozen and inside the
            % picture
            if((i>0)&&(j>0)&&(i<=size(F,1))&&(j<=size(F,2))&&(Frozen(i,j)==0))
                Tt=1/max(F(i,j),eps);
                Ty=1;
                % Update distance in neigbour list or add to neigbour list
                if(T(i,j)>0)
                    if(neg_list(1,T(i,j))>Tt)
                        neg_list(1,T(i,j))=Tt;
                    end
                    if(Ed)
                        neg_list(4,T(i,j))=min(Ty,neg_list(4,T(i,j)));
                    end
                else
                    neg_pos=neg_pos+1;
                    % If running out of memory at a new block
                    if(neg_pos>neg_free), neg_free = neg_free +100000; neg_list(1,neg_free)=0; end
                    if(Ed)
                        neg_list(:,neg_pos)=[Tt;i;j;Ty];
                    else
                        neg_list(:,neg_pos)=[Tt;i;j];
                    end
                    T(i,j)=neg_pos;
                end
            end
        end
    end

    % Loop through all pixels of the image
    for itt=1:numel(F)
        % Get the pixel from narrow list (boundary list) with smallest
        % distance value and set it to current pixel location
        [t,index]=min(neg_list(1,1:neg_pos));
        if(neg_pos==0), break; end
        x=neg_list(2,index); y=neg_list(3,index);
        Frozen(x,y)=1;
        T(x,y)=neg_list(1,index);

        if(Ed), Y(x,y)=neg_list(4,index); end

        % Remove min value by replacing it with the last value in the array
        if(index<neg_pos),
            neg_list(:,index)=neg_list(:,neg_pos);
            x2=neg_list(2,index); y2=neg_list(3,index);
            T(x2,y2)=index; 
        end
        neg_pos =neg_pos-1;

        % Loop through all 4 neighbours of current pixel
        for k=1:4,
            % Location of neighbour
            i=x+ne(k,1); j=y+ne(k,2);

            % Check if current neighbour is not yet frozen and inside the
            % picture
            if((i>0)&&(j>0)&&(i<=size(F,1))&&(j<=size(F,2))&&(Frozen(i,j)==0))

                Tt=CalculateDistance(T,F(i,j),size(F),i,j,usesecond,usecross,Frozen);
                if(Ed)
                    Ty=CalculateDistance(Y,1,size(F),i,j,usesecond,usecross,Frozen);
                end

                % Update distance in neigbour list or add to neigbour list
                if(T(i,j)>0)
                    neg_list(1,T(i,j))=min(Tt,neg_list(1,T(i,j)));
                    if(Ed)
                        neg_list(4,T(i,j))=min(Ty,neg_list(4,T(i,j)));
                    end
                else
                    neg_pos=neg_pos+1;
                    % If running out of memory at a new block
                    if(neg_pos>neg_free), neg_free = neg_free +100000; neg_list(1,neg_free)=0; end
                    if(Ed)
                        neg_list(:,neg_pos)=[Tt;i;j;Ty];
                    else
                        neg_list(:,neg_pos)=[Tt;i;j];
                    end
                    T(i,j)=neg_pos;
                end
            end
        end
    end
end

function Tt=CalculateDistance(T,Fij,sizeF,i,j,usesecond,usecross,Frozen)
    % Boundary and frozen check -> current patch
    Tpatch=inf(5,5);
    for nx=-2:2
        for ny=-2:2
            in=i+nx; jn=j+ny;
            if((in>0)&&(jn>0)&&(in<=sizeF(1))&&(jn<=sizeF(2))&&(Frozen(in,jn)==1))
                Tpatch(nx+3,ny+3)=T(in,jn);
            end
        end
    end

    % The values in order is 0 if no neighbours in that direction
    % 1 if 1e order derivatives is used and 2 if second order
    % derivatives are used
    Order=zeros(1,4);

    % Make 1e order derivatives in x and y direction
    Tm(1) = min( Tpatch(2,3) , Tpatch(4,3)); if(isfinite(Tm(1))), Order(1)=1; end
    Tm(2) = min( Tpatch(3,2) , Tpatch(3,4)); if(isfinite(Tm(2))), Order(2)=1; end
    % Make 1e order derivatives in cross directions
    if(usecross)
        Tm(3) = min( Tpatch(2,2) , Tpatch(4,4)); if(isfinite(Tm(3))), Order(3)=1; end
        Tm(4) = min( Tpatch(2,4) , Tpatch(4,2)); if(isfinite(Tm(4))), Order(4)=1; end
    end

    % Make 2e order derivatives
    if(usesecond)
        Tm2=zeros(1,4);
        % pixels with a pixeldistance 2 from the center must be
        % lower in value otherwise use other side or first order
        ch1=(Tpatch(1,3)<Tpatch(2,3))&&isfinite(Tpatch(2,3)); ch2=(Tpatch(5,3)<Tpatch(4,3))&&isfinite(Tpatch(4,3));

        if(ch1&&ch2),Tm2(1) =min( (4*Tpatch(2,3)-Tpatch(1,3))/3 , (4*Tpatch(4,3)-Tpatch(5,3))/3);  Order(1)=2;
        elseif(ch1), Tm2(1) =(4*Tpatch(2,3)-Tpatch(1,3))/3; Order(1)=2;
        elseif(ch2), Tm2(1) =(4*Tpatch(4,3)-Tpatch(5,3))/3; Order(1)=2;
        end

        ch1=(Tpatch(3,1)<Tpatch(3,2))&&isfinite(Tpatch(3,2)); ch2=(Tpatch(3,5)<Tpatch(3,4))&&isfinite(Tpatch(3,4));

        if(ch1&&ch2),Tm2(2) =min( (4*Tpatch(3,2)-Tpatch(3,1))/3 , (4*Tpatch(3,4)-Tpatch(3,5))/3); Order(2)=2;
        elseif(ch1), Tm2(2)=(4*Tpatch(3,2)-Tpatch(3,1))/3; Order(2)=2;
        elseif(ch2), Tm2(2)=(4*Tpatch(3,4)-Tpatch(3,5))/3; Order(2)=2;
        end

        if(usecross)
            ch1=(Tpatch(1,1)<Tpatch(2,2))&&isfinite(Tpatch(2,2)); ch2=(Tpatch(5,5)<Tpatch(4,4))&&isfinite(Tpatch(4,4));
            if(ch1&&ch2),Tm2(3) =min( (4*Tpatch(2,2)-Tpatch(1,1))/3 , (4*Tpatch(4,4)-Tpatch(5,5))/3); Order(3)=2;
            elseif(ch1), Tm2(3)=(4*Tpatch(2,2)-Tpatch(1,1))/3; Order(3)=2;
            elseif(ch2), Tm2(3)=(4*Tpatch(4,4)-Tpatch(5,5))/3; Order(3)=2;
            end

            ch1=(Tpatch(1,5)<Tpatch(2,4))&&isfinite(Tpatch(2,4)); ch2=(Tpatch(5,1)<Tpatch(4,2))&&isfinite(Tpatch(4,2));
            if(ch1&&ch2),Tm2(4) =min( (4*Tpatch(2,4)-Tpatch(1,5))/3 , (4*Tpatch(4,2)-Tpatch(5,1))/3); Order(4)=2;
            elseif(ch1), Tm2(4)=(4*Tpatch(2,4)-Tpatch(1,5))/3; Order(4)=2;
            elseif(ch2), Tm2(4)=(4*Tpatch(4,2)-Tpatch(5,1))/3; Order(4)=2;
            end
        end
    else
        Tm2=zeros(1,4);
    end

    % Calculate the distance using x and y direction
    Coeff = [0 0 -1/(max(Fij^2,eps))];
    for t=1:2;
        switch(Order(t))
            case 1,
                Coeff=Coeff+[1 -2*Tm(t) Tm(t)^2];
            case 2,
                Coeff=Coeff+[1 -2*Tm2(t) Tm2(t)^2]*(2.2500);
        end
    end

    Tt=roots(Coeff); Tt=max(Tt);
    % Calculate the distance using the cross directions
    if(usecross)
        Coeff = Coeff + [0 0 -1/(max(Fij^2,eps))];
        for t=3:4;
            switch(Order(t))
                case 1,
                    Coeff=Coeff+0.5*[1 -2*Tm(t) Tm(t)^2];
                case 2,
                    Coeff=Coeff+0.5*[1 -2*Tm2(t) Tm2(t)^2]*(2.2500);
            end
        end
        Tt2=roots(Coeff); Tt2=max(Tt2);
        % Select minimum distance value of both stensils
        if(~isempty(Tt2)), Tt=min(Tt,Tt2); end
    end

    % Upwind condition check, current distance must be larger
    % then direct neighbours used in solution 
    DirectNeigbInSol=Tm(isfinite(Tm));
    if(nnz(DirectNeigbInSol>=Tt)>0) % Will this ever happen?
        Tt=min(DirectNeigbInSol)+(1/(max(Fij,eps)));
    end
end

function z=roots(Coeff)
    a=Coeff(1); b=Coeff(2); c=Coeff(3); d=max((b*b)-4.0*a*c,0);
    if(a~=0)
        z(1)= (-b - sqrt(d)) / (2.0*a);
        z(2)= (-b + sqrt(d)) / (2.0*a);
    else 
        z(1)= (2.0*c)/(-b - sqrt(d));
        z(2)= (2.0*c)/(-b + sqrt(d));
    end
end

function path = compute_path(arrival_time_map, end_points, ray_grid_step_size)

    % define the output to be empty cell structure.
    path = {};

    % get size of end points. these points is ray end points but path start points.
    [end_point_dimensions, end_points_num] = size(end_points);

    % check input parameters.
    if end_point_dimensions ~=2
        display('input error!!!');
        return;
    end


    % rewrite inf in arrival_time_map to its maximum value
    I = find(isinf(arrival_time_map));
    J = find(~isinf(arrival_time_map));
    A1 = arrival_time_map;
    A1(I) = max(arrival_time_map(J));

    % calculate gradient
    grad = zeros(256, 256,2);
    [grad(:,:,1),grad(:,:,2)] = gradient(A1);

    n = sum(grad.^2, 3);
    n(find(n<eps)) = 1;
    grad = -grad ./ repmat(n,[1,1,2]);


    path = stream2(grad(:,:,1),grad(:,:,2),end_points(2,:),end_points(1,:), [ray_grid_step_size, 1000]);  % last parameter are ray grid step size (unit is pixel), and largest step number.
    for point_i = 1:length(path)
        path{point_i} = (path{point_i}(:,2:-1:1))';


        d = sum ((path{point_i} - repmat(path{point_i}(:,end),[1,length(path{point_i})])).^2);
        T = max(d)/300^2;
        I = find(d<T);
        if  not(isempty(I))
            path{point_i} = path{point_i}(:,1:I(1));
        end



    end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    coefficient projection matrix calculation for Curve (Bent) ray reconstruction method using virtual receiver technique. 
%
%   Input:
%       ring_diameter       transducer ring diameter in meter.  
%       element_num         transducer elements number.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
%       curr_image          current image is the base for ray tracing.
%       pjt_data            travel time difference map. 
%                           1st D. is receiver, 
%                           2nd D. is emitter, 
%       ray_grid_step_size  grid step size for ray tracing(unit is meter).
% 
%   Output:
%       a_ray_matrix:       ray path weight matrix (coefficient projection matrix) without hanning window.
%       w_ray_matrix:       same to last output but with hanning window.
%       pjt_data_vr:        travel time difference map for vitual recier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_ray_matrix, w_ray_matrix, pjt_data_vr] = curve_ray_tracing_CVR(ring_diameter, element_num, image_size_p, image_size_m, curr_image, pjt_data, ray_grid_step_size, speed_in_water)
    
    % ray grid step size
    ds = ray_grid_step_size;

    % calculate refractive index map using inputed curr_image with sound speed information.
    n = speed_in_water./curr_image;                % refractive index map
    h = fspecial('average', [3 3]);     % average filter kernel
    n = filter2(h,n);                   % average smoothed refractive index map.

    % calculate each element position, angle of start position is 0;
    theta_rg = linspace(0, 2*pi*(element_num-1)/element_num, element_num);
    ring_radius = ring_diameter/2;                      % radius of the ring.        
    element_pos_x = ring_radius * cos(theta_rg);        % element position x
    element_pos_y = ring_radius * sin(theta_rg);        % elmeent position y.
    
    % extent elements anlge for angles interpolation.
    theta_rg_extent = [theta_rg-2*pi, theta_rg, theta_rg+2*pi];
    
    % allocate space to save traced ray and modified travle time difference map.
    ray_grid_pos_x = ones(3000,element_num,element_num).*inf;   % saving all ray grid steps position of each traced curve ray.ˆäãC³1000¨3000
    ray_grid_pos_y = ray_grid_pos_x;
    pjt_data_vr = zeros(size(pjt_data));                        % saving modified projection data. 
    
    
    % for loop for each transmit.
    for transmit_i = 1:element_num
        
        % inital_shoot_angle for current transmit.
        shoot_angle = angle(complex(element_pos_x-element_pos_x(transmit_i), element_pos_y-element_pos_y(transmit_i)));

        % allocate space to record arrival postion of each curve ray.
        arrival_pos_x = zeros(1,element_num);
        arrival_pos_y = zeros(1,element_num);
        
        
        % loop for each receiver or for each shoot angle.
        for shoot_angle_i = 1:element_num
            
            % when receiver and transmit are same elements
            if shoot_angle_i == transmit_i
                
                % arrival position will be current transmit or receiver
                arrival_pos_x(shoot_angle_i) = element_pos_x(transmit_i);
                arrival_pos_y(shoot_angle_i) = element_pos_y(transmit_i);
                
                % first poition for ray tracing which is transmitter element position.
                ray_grid_pos_x(1,shoot_angle_i,transmit_i) = element_pos_x(transmit_i);
                ray_grid_pos_y(1,shoot_angle_i,transmit_i) = element_pos_y(transmit_i);
                
                continue;
            end
            
           
            theta1 = shoot_angle(shoot_angle_i);    % current shoot angle.
            
            % current shoot angle tracing.
            [single_ray_grid_pos_x, single_ray_grid_pos_y, single_ray_tof] = single_curve_ray_tracing([element_pos_x(transmit_i), element_pos_y(transmit_i)], theta1, ring_radius, image_size_p, image_size_m, n, ds,speed_in_water);
            
            % record arrival position of current ray
            arrival_pos_x(shoot_angle_i) = single_ray_grid_pos_x(end);
            arrival_pos_y(shoot_angle_i) = single_ray_grid_pos_y(end);
            
            % record each grid steps of current ray
            ray_grid_pos_x(1:length(single_ray_grid_pos_x),shoot_angle_i,transmit_i) = single_ray_grid_pos_x';
            ray_grid_pos_y(1:length(single_ray_grid_pos_y),shoot_angle_i,transmit_i) = single_ray_grid_pos_y';
        end
        
        % angles of curve ray arrival positions in polar coordinates 
        arrival_angle = angle(complex(arrival_pos_x,arrival_pos_y));
        arrival_angle(find(arrival_angle<0)) = arrival_angle(find(arrival_angle<0)) + 2*pi;
        
        % measured tofs of virtual receivers obtained by interpolation.
        pjt_data_vr(:,transmit_i) = interp1(theta_rg_extent, repmat(pjt_data(:,transmit_i),[3,1]), arrival_angle, 'linear');
        
    end

    % get coefficient projection matrix.
    [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m);
    
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    coefficient projection matrix calculation for Curve (Bent) ray reconstruction method using Newton-Raphson method for ray linking. 
%
%   Input:
%       ring_diameter       transducer ring diameter in meter.  
%       element_num         transducer elements number.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
%       curr_image          current image is the base for ray tracing.
%       ray_grid_step_size  grid step size for ray tracing(unit is meter).
% 
%   Output:
%       a_ray_matrix:       ray path weight matrix (coefficient projection matrix) without hanning window.
%       w_ray_matrix:       same to last output but with hanning window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_ray_matrix, w_ray_matrix] = curve_ray_tracing_CNR(ring_diameter, element_num, image_size_p, image_size_m, curr_image, ray_grid_step_size, speed_in_water)
    
    % error tolerance parameter.
    epsiron = 1.e-4;                                        % 0.1 mm error tollerance.
    error_tol_angle = epsiron/(pi*ring_diameter) *2*pi;     % same error tollerance in arrival angles.
    
    
    % curve ray tracing step size.
    ds = ray_grid_step_size;
    
    % calculate refractive index map using inputed curr_image with sound speed information.
    n = speed_in_water./curr_image;                % refractive index map
    h = fspecial('average', [3 3]);     % average filter kernel
    n = filter2(h,n);                   % average smoothed refractive index map.

    % calculate each element position, angle of start position is 0*pi;
    theta_rg = linspace(0, 2*pi*(element_num-1)/element_num , element_num);     % element angles in the ring.
    ring_radius = ring_diameter/2;                                              % radius of the ring.        
    element_pos_x = ring_radius * cos(theta_rg);                                % element position x
    element_pos_y = ring_radius * sin(theta_rg);                                % elmeent position y
    
    % allocate space for saveing each traced curve ray. Note: maximum ray lenght is 1000 ray steps.
    ray_grid_pos_x = ones(1000,element_num,element_num).*inf;   % x position
    ray_grid_pos_y = ray_grid_pos_x;                            % y position

    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%     counter_failed_ray_num = 0;
%     counter_traced_ray_num = 0;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
    
    
    %loop for each transmit element
    for transmit_i = 1:element_num

        % inital_shoot_angle for current transmit.
        initial_shoot_angle = linspace(theta_rg(transmit_i)+pi/2, theta_rg(transmit_i)+pi*3/2-pi/element_num, element_num);
        initial_shoot_angle = circshift(initial_shoot_angle,[0,transmit_i-1]);              % shift for coordinate order.
        
        % loop for each receive element
        for receive_i = 1:element_num

            % when both transmiter and receiver are same elements.
            if transmit_i == receive_i
                continue;
            end
            
            
            % prepare iteration parameters.
            iteration_max_num = 20;                                                 % maximum iteration number for find a single curve ray
            iteration_arrival_pos_angle = ones(1,iteration_max_num)*inf;            % arrival angle position for each iteration steps
            iteration_arrival_pos_angle_abs_errors = iteration_arrival_pos_angle;   % arrival angle absolute errors.
            iteration_shoot_angle = iteration_arrival_pos_angle;                    % shoot angle for each iteration steps.
            
            
            % ring linking problem should be solved in following iteration
            % the function should be "arrival_angle = a * shoot_angle + b", a is slop, it is written by k in code.
            % thus try to solve function "a * shoot_angle + b - arrival_angle = 0", which is "f(x) = 0" shoot_angle is x.
            % "f(x)=0" can by solved by Newton Raphson iteration
            % next loop is Newton Raphson iteration.
            for iteration_i = 1:iteration_max_num
            
                % in first iteration, initial shoot angle is employed
                if iteration_i ==1
                    theta1 = initial_shoot_angle(receive_i);    % theta1 is current shoot angle. 
                
                % not first iteration
                else
                    
                    % in 2nd iteration, slop can not be calculated.   
                    if iteration_i == 2
                        k = 2;      % slop value, when homogenous.
                        
                    % not 1st or 2nd iteration, slop can be calculated bu prevous two point.    
                    else
                        
                        % arrival angle different of prevous two points
                        arrival_angle_diff = iteration_arrival_pos_angle(iteration_i-2) - iteration_arrival_pos_angle(iteration_i-1);
                        if abs(arrival_angle_diff) > 2*pi - abs(arrival_angle_diff)  
                            arrival_angle_diff = angle(complex(cos(iteration_arrival_pos_angle(iteration_i-2)), sin(iteration_arrival_pos_angle(iteration_i-2)))) - angle(complex(cos(iteration_arrival_pos_angle(iteration_i-1)), sin(iteration_arrival_pos_angle(iteration_i-1))));    
                        end
                        
                        % k is very important, we may have to remove, very small k, since it may cause iterations diverge.
                        k = arrival_angle_diff / (iteration_shoot_angle(iteration_i-2) - iteration_shoot_angle(iteration_i-1));
                        
                    end
                    
                    % arrival angle errors
                    arrival_angle_error = arrival_theta-theta_rg(receive_i);
                    if abs(arrival_angle_error) > 2*pi - abs(arrival_angle_error)
                        arrival_angle_error = angle(complex(cos(arrival_theta), sin(arrival_theta))) - angle(complex(cos(theta_rg(receive_i)), sin(theta_rg(receive_i))));   
                    end
                    
                    theta1 = theta1 - 0.8 * arrival_angle_error/k;  % Newton Raphson iteration. 0.8 is for relaxation.

                end

                % trace curve ray with current transmit and shoot angle. 
                [single_ray_grid_pos_x, single_ray_grid_pos_y, single_ray_tof] = single_curve_ray_tracing([element_pos_x(transmit_i), element_pos_y(transmit_i)], theta1, ring_radius, image_size_p, image_size_m, n, ds,speed_in_water);

                % arrival position angle in polar coordinate system.
                arrival_theta = angle(complex(single_ray_grid_pos_x(end),single_ray_grid_pos_y(end)));
                if arrival_theta < 0
                    arrival_theta = arrival_theta + 2*pi;
                end
                
                % record current shoot angle and arrival position angle.
                iteration_arrival_pos_angle(iteration_i) = arrival_theta;
                iteration_shoot_angle(iteration_i) = theta1;

                % if linked keep path.
                iteration_arrival_pos_angle_abs_errors(iteration_i) = abs(arrival_theta - theta_rg(receive_i));   % absolute errors of iteration arrival position angle 
                iteration_arrival_pos_angle_abs_errors(iteration_i) = min(iteration_arrival_pos_angle_abs_errors(iteration_i), 2*pi-iteration_arrival_pos_angle_abs_errors(iteration_i));
                
                % current traced curve ray arrived expected position.
                if iteration_arrival_pos_angle_abs_errors(iteration_i) < error_tol_angle        
                    
                    % record each step position of traced curve ray.
                    ray_grid_pos_x(1:length(single_ray_grid_pos_x),receive_i,transmit_i) = single_ray_grid_pos_x';    
                    ray_grid_pos_y(1:length(single_ray_grid_pos_y),receive_i,transmit_i) = single_ray_grid_pos_y'; 
                    break; 
                
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record number of ray tracing
%                     counter_traced_ray_num = counter_traced_ray_num+iteration_i;
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
                    
                % curve ray failed to arrive expeted position untile maximum iteration step.
                elseif iteration_i == iteration_max_num
                                       
                    % straight ray tracing since curve ray tracing failed.
                    N = sqrt((element_pos_x(receive_i)-element_pos_x(transmit_i))^2+(element_pos_y(receive_i)-element_pos_y(transmit_i))^2)/ds;
                    pos_x = linspace(element_pos_x(transmit_i),element_pos_x(receive_i),N);
                    pos_y = linspace(element_pos_y(transmit_i),element_pos_y(receive_i),N);
                    ray_grid_pos_x(1:length(pos_x),receive_i,transmit_i) = pos_x';
                    ray_grid_pos_y(1:length(pos_y),receive_i,transmit_i) = pos_y';

%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%                     counter_failed_ray_num = counter_failed_ray_num+1;
%                     counter_traced_ray_num = counter_traced_ray_num+iteration_i;
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
                    
                end
            end 
        end
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%     display(['traced ray number is ',num2str(counter_traced_ray_num),' and failed: ', num2str(counter_failed_ray_num)]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
    
    % get Coefficient matrix by curve ray tracing result.
    [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    coefficient projection matrix calculation for Curve (Bent) ray reconstruction method using Bisection method for ray linking. 
%
%   Input:
%       ring_diameter       transducer ring diameter in meter.  
%       element_num         transducer elements number.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
%       curr_image          current image is the base for ray tracing.
%       ray_grid_step_size  grid step size for ray tracing(unit is meter).
% 
%   Output:
%       a_ray_matrix:       ray path weight matrix (coefficient projection matrix) without hanning window.
%       w_ray_matrix:       same to last output but with hanning window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_ray_matrix, w_ray_matrix] = curve_ray_tracing_CBS(ring_diameter, element_num, image_size_p, image_size_m, curr_image, ray_grid_step_size,speed_in_water)
    
    % error tolerance parameter.
    epsiron = 1.e-4;                                        % 0.1 mm error tollerance.
    error_tol_angle = epsiron/(pi*ring_diameter) *2*pi;     % same error tollerance in arrival angles.
    
    
    % curve ray tracing step size.
    ds = ray_grid_step_size;   % ray grid size, generally should be half of recon image resolution.

    % calculate refractive index map using inputed curr_image with sound speed information.
    n = speed_in_water./curr_image;                % refractive index map
    h = fspecial('average', [3 3]);     % average filter kernel
    n = filter2(h,n);                   % average smoothed refractive index map.

    % calculate each element position, start position is 0*pi;
    theta_rg = linspace(0, 2*pi*(element_num-1)/element_num , element_num);     % element angles in the ring.
    ring_radius = ring_diameter/2;                                              % radius of the ring.        
    element_pos_x = ring_radius * cos(theta_rg);                                % element position x
    element_pos_y = ring_radius * sin(theta_rg);                                % elmeent position y
    
    % allocate space for saveing each traced curve ray. Note: maximum ray lenght is 1000 ray steps.
    ray_grid_pos_x = ones(1000,element_num,element_num).*inf;   % x position
    ray_grid_pos_y = ray_grid_pos_x;                            % y position
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%     counter_failed_ray_num = 0;
%     counter_traced_ray_num = 0;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
    
    
    %loop for each transmit element
    for transmit_i = 1:element_num

        % inital_shoot_angle for current transmit.
        initial_shoot_angle = linspace(theta_rg(transmit_i)+pi/2, theta_rg(transmit_i)+pi*3/2, element_num+1);
        initial_shoot_angle(1)   = initial_shoot_angle(1)   + 0.1*pi/element_num;
        initial_shoot_angle(end) = initial_shoot_angle(end) - 0.1*pi/element_num;

        % allocate space to save angles of arrival position
        initial_arrival_pos_angle = ones(1,numel(initial_shoot_angle)).*inf;
        
        % loop for each initial shoot angle
        for initial_shoot_angle_i = 1:numel(initial_shoot_angle)
            
            % current initial shoot angle
            theta1 = initial_shoot_angle(initial_shoot_angle_i);
            
            % trace curve ray with current transmit and shoot angle. 
            [single_ray_grid_pos_x, single_ray_grid_pos_y, single_ray_tof] = single_curve_ray_tracing([element_pos_x(transmit_i), element_pos_y(transmit_i)], theta1, ring_radius, image_size_p, image_size_m, n, ds,speed_in_water);
            
            % angle of arrival position
            arrival_theta = angle(complex(single_ray_grid_pos_x(end),single_ray_grid_pos_y(end)));
            if arrival_theta < 0
                arrival_theta = arrival_theta + 2*pi;
            end
            
            % record current shoot angle and arrival position angle.
            initial_arrival_pos_angle(initial_shoot_angle_i) = arrival_theta;            
        end
        
        
        
        % sort angle of arrival position for all initial shoot angles.
        [initial_sorted_arrival_angle, initial_sorted_index] = sort(initial_arrival_pos_angle);
        
        % loop for each virtual receive element
        for receive_i = 1:element_num
            
            % when both transmiter and receiver are same elements.
            if transmit_i == receive_i
                continue;
            end
            
            % prepare iteration parameters.
            iteration_max_num = 20;                                                 % maximum iteration number for find a single curve ray     
            iteration_arrival_pos_angle = ones(1,iteration_max_num)*inf;            % arrival angle position for each iteration steps
            iteration_shoot_angle = iteration_arrival_pos_angle;                    % shoot angle for each iteration steps.

            % select two start point from angles of arrival position from initial shoot angles.
            [max_value,max_index] = max(initial_sorted_arrival_angle >theta_rg(receive_i));

            % when there are arrivale angles one larger than required and one less than required. both angles will be initial starting of bisection mtehod.  
            if any(initial_arrival_pos_angle>theta_rg(receive_i)) && any(initial_arrival_pos_angle<theta_rg(receive_i))
                iteration_shoot_angle(1) = initial_shoot_angle(initial_sorted_index(max_index));
                iteration_shoot_angle(2) = initial_shoot_angle(initial_sorted_index(max_index-1));

                iteration_arrival_pos_angle(1) = initial_arrival_pos_angle(initial_sorted_index(max_index));
                iteration_arrival_pos_angle(2) = initial_arrival_pos_angle(initial_sorted_index(max_index-1));

            % if there is no goog initial starting.
            else
                iteration_shoot_angle(1) = initial_shoot_angle(initial_sorted_index(end));
                iteration_shoot_angle(2) = initial_shoot_angle(initial_sorted_index(1));

                iteration_arrival_pos_angle(1) = initial_arrival_pos_angle(initial_sorted_index(end));
                iteration_arrival_pos_angle(2) = initial_arrival_pos_angle(initial_sorted_index(1));
            end
            
            % record active shoot angle and their arrival angles. both active shoot angle will be used to calculate next new shoot angle using bisection methdo.
            active_mark = [1,2];
            
            % iteration using Bisection method 
            for iteration_i = 3:iteration_max_num
            
                % in third iteration, there are just two prevous points, thus they are used to calculation next new shoot angle. 
                if iteration_i == 3
                    current_shoot_angle = angle(complex(cos(iteration_shoot_angle(1))+cos(iteration_shoot_angle(2)), sin(iteration_shoot_angle(1))+sin(iteration_shoot_angle(2))));
                else
                    
                    % when their are up and low arrivale angles, they will be used to calculate next new shoot angle.
                    if any(iteration_arrival_pos_angle(active_mark)>theta_rg(receive_i)) && any(iteration_arrival_pos_angle(active_mark)<theta_rg(receive_i))
                        
                        % when last new shoot anvle's arrival angle larger than required arrival angle.
                        if iteration_arrival_pos_angle(iteration_i-1) > theta_rg(receive_i)
                            
                            % change active mark.
                            if iteration_arrival_pos_angle(active_mark(1)) > theta_rg(receive_i)
                                active_mark(1) = iteration_i-1;
                            else
                                active_mark(2) = iteration_i-1;
                            end
                        
                        % when small
                        else
                            
                            % change active mark.
                            if iteration_arrival_pos_angle(active_mark(1)) < theta_rg(receive_i)
                                active_mark(1) = iteration_i-1;
                            else
                                active_mark(2) = iteration_i-1;
                            end
                        end
                        
                        % claulcate next new shoot angle
                        current_shoot_angle = angle(complex(cos(iteration_shoot_angle(active_mark(1)))+cos(iteration_shoot_angle(active_mark(2))), sin(iteration_shoot_angle(active_mark(1)))+sin(iteration_shoot_angle(active_mark(2)))));
                    
                    % when there is no both larger and less arrival angle, the arrival angle should locate very close to 0 degree.
                    else
                        
                        % do not use arrival angle, use sin of arrival angle.
                        if sin(iteration_arrival_pos_angle(iteration_i-1)) > sin(theta_rg(receive_i))
                        
                            if sin(iteration_arrival_pos_angle(active_mark(1))) > sin(theta_rg(receive_i))
                                active_mark(1) = iteration_i-1;
                            else
                                active_mark(2) = iteration_i-1;
                            end
                        
                        % when last new arrival angle less than required arrival angle. 
                        else
                            if sin(iteration_arrival_pos_angle(active_mark(1))) < sin(theta_rg(receive_i))
                                active_mark(1) = iteration_i-1;
                            else
                                active_mark(2) = iteration_i-1;
                            end
                            
                        end
                        
                        % calucate new shoot angle.
                        current_shoot_angle = angle(complex(cos(iteration_shoot_angle(active_mark(1)))+cos(iteration_shoot_angle(active_mark(2))), sin(iteration_shoot_angle(active_mark(1)))+sin(iteration_shoot_angle(active_mark(2)))));
                        
                    end
                    
                end
                

                % trace curve ray with current transmit and shoot angle. 
                [single_ray_grid_pos_x, single_ray_grid_pos_y, single_ray_tof] = single_curve_ray_tracing([element_pos_x(transmit_i), element_pos_y(transmit_i)], current_shoot_angle, ring_radius, image_size_p, image_size_m, n, ds,speed_in_water);

                
                % arrival position angle in polar coordinate system.
                current_arrival_angle = angle(complex(single_ray_grid_pos_x(end),single_ray_grid_pos_y(end)));
                if current_arrival_angle < 0
                    current_arrival_angle = current_arrival_angle + 2*pi;
                end
               
                iteration_shoot_angle(iteration_i) = current_shoot_angle;
                iteration_arrival_pos_angle(iteration_i) = current_arrival_angle;
                
                
                % current traced curve ray arrived expected position.
                if min(abs(current_arrival_angle - theta_rg(receive_i)), 2*pi-abs(current_arrival_angle - theta_rg(receive_i))) < error_tol_angle        
                    
                    % record each step position of traced curve ray.
                    ray_grid_pos_x(1:length(single_ray_grid_pos_x),receive_i,transmit_i) = single_ray_grid_pos_x';    
                    ray_grid_pos_y(1:length(single_ray_grid_pos_y),receive_i,transmit_i) = single_ray_grid_pos_y'; 
                     
                    
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%                     counter_traced_ray_num = counter_traced_ray_num+iteration_i;
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
                    
                    break;
                    
                % curve ray failed to arrive expeted position untile maximum iteration step.
                elseif iteration_i == iteration_max_num
          
                    
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%                     counter_failed_ray_num = counter_failed_ray_num+1;
%                     counter_traced_ray_num = counter_traced_ray_num+iteration_i;
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing
                    
                    
                    % straight ray tracing since curve ray tracing failed.
                    N = sqrt((element_pos_x(receive_i)-element_pos_x(transmit_i))^2+(element_pos_y(receive_i)-element_pos_y(transmit_i))^2)/ds;
                    pos_x = linspace(element_pos_x(transmit_i),element_pos_x(receive_i),N);
                    pos_y = linspace(element_pos_y(transmit_i),element_pos_y(receive_i),N);
                    ray_grid_pos_x(1:length(pos_x),receive_i,transmit_i) = pos_x';
                    ray_grid_pos_y(1:length(pos_y),receive_i,transmit_i) = pos_y';
                end
            end 
        end
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for record numbers of ray tracing and failed ray linking, respectively.
%     display(['traced ray number is ',num2str(counter_traced_ray_num),' and failed: ', num2str(counter_failed_ray_num)]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%testing
    
    % get Coefficient matrix by curve ray tracing result.
    [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m);
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    coefficient projection matrix calculation for Curve (Bent) ray reconstruction method using burte force method for ray linking. 
%
%   Input:
%       ring_diameter       transducer ring diameter in meter.  
%       element_num         transducer elements number.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
%       curr_image          current image is the base for ray tracing.
%       ray_grid_step_size  grid step size for ray tracing(unit is meter).
% 
%   Output:
%       a_ray_matrix:       ray path weight matrix (coefficient projection matrix) without hanning window.
%       w_ray_matrix:       same to last output but with hanning window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_ray_matrix, w_ray_matrix] = curve_ray_tracing_CBF(ring_diameter, element_num, image_size_p, image_size_m, curr_image,ray_grid_step_size, speed_in_water)
    
    % error tolerance parameter.
    epsiron = 1.e-4;                                        % 0.1 mm error tollerance.
    %error_tol_angle = epsiron/(pi*ring_diameter) *2*pi;     % same error tollerance in arrival angles.
    
    
    % curve ray tracing step size.
    %ds = 1540/1.e6/4;                  % 0.385mm ray tracing step size.
    ds = ray_grid_step_size;   % ray grid size, generally should be half of recon image resolution.

    % calculate refractive index map using inputed curr_image with sound speed information.
    n = speed_in_water./curr_image;                % refractive index map
    h = fspecial('average', [3 3]);     % average filter kernel
    n = filter2(h,n);                   % average smoothed refractive index map.

    % calculate each element position, start position is 0*pi;
    theta_rg = linspace(0, 2*pi*(element_num-1)/element_num , element_num);     % element angles in the ring.
    ring_radius = ring_diameter/2;                                              % radius of the ring.        
    element_pos_x = ring_radius * cos(theta_rg);                                % element position x
    element_pos_y = ring_radius * sin(theta_rg);                                % elmeent position y
    
    % allocate space for saveing each traced curve ray. Note: maximum ray lenght is 1000 ray steps.
    ray_grid_pos_x = ones(1000,element_num,element_num).*inf;   % x position
    ray_grid_pos_y = ray_grid_pos_x;                            % y position

    
    %loop for each transmit element
    for transmit_i = 1:element_num

        % inital_shoot_angle for current transmit.
        shoot_angle = linspace(theta_rg(transmit_i)+pi/2, theta_rg(transmit_i)+pi*3/2-pi/element_num, element_num*2^9);
        
        % allocate space for saving shortest tof of each receiver in muti ray path.
        ch_shoot_tof = zeros(1,element_num);
        for shoot_angle_i = 1 : length(shoot_angle)
       
            % trace curve ray with current transmit and shoot angle. 
            [single_ray_grid_pos_x, single_ray_grid_pos_y, single_ray_tof] = single_curve_ray_tracing([element_pos_x(transmit_i), element_pos_y(transmit_i)], shoot_angle(shoot_angle_i), ring_radius, image_size_p, image_size_m, n, ds,speed_in_water);
            
            % compare arrival position with element position
            arrival_pos = repmat([single_ray_grid_pos_x(end) ; single_ray_grid_pos_y(end)] , 1 , element_num);  % arrival position of current curve ray tracing.
            arrival_pos_diff = arrival_pos - [element_pos_x; element_pos_y];                                    % difference between arrival position and each receiver position.
            arrival_pos_error = sqrt(arrival_pos_diff(1,:).^2+(arrival_pos_diff(2,:).^2));                      % absolute error between arrival position and each receiver position.
            
            [min_value, min_index] = min(arrival_pos_error);        %  the closest receiver with its index.
            
            % if current curve ray arrived the closest receiver.
            if min_value < epsiron

                if ch_shoot_tof(min_index)==0                   % first time ray arrived this receiver.
                    ray_grid_pos_x(1:length(single_ray_grid_pos_x),min_index,transmit_i) = single_ray_grid_pos_x';      % record each ray steps.
                    ray_grid_pos_y(1:length(single_ray_grid_pos_y),min_index,transmit_i) = single_ray_grid_pos_y';                        
                    ch_shoot_tof(min_index) = single_ray_tof;                                                           % record ray tof

                elseif ch_shoot_tof(min_index)>single_ray_tof   % not first time
                    ray_grid_pos_x(:,min_index,transmit_i) = inf;                                                       % record each ray steps.
                    ray_grid_pos_x(1:length(single_ray_grid_pos_x),min_index,transmit_i) = single_ray_grid_pos_x';
                    ray_grid_pos_y(:,min_index,transmit_i) = inf;
                    ray_grid_pos_y(1:length(single_ray_grid_pos_y),min_index,transmit_i) = single_ray_grid_pos_y';                        
                    ch_shoot_tof(min_index) = single_ray_tof;                                                           % record ray tof.
                end             
            end      
        end
        
        mis_index = find(ch_shoot_tof == 0);        % find recevers which has no curve ray arrived.
        mis_num = numel(mis_index);
        
        if  mis_num > 0
            for straight_ray_i = 1:mis_num

                % if this receiver is not same element to transmit
                if mis_index(straight_ray_i) ~= transmit_i

                    % saving straight line for curent pair of transmiter and receiver.
                    N = sqrt((element_pos_x(mis_index(straight_ray_i))-element_pos_x(transmit_i))^2+(element_pos_y(mis_index(straight_ray_i))-element_pos_y(transmit_i))^2)/ds;
                    pos_x = linspace(element_pos_x(transmit_i),element_pos_x(mis_index(straight_ray_i)),N);
                    pos_y = linspace(element_pos_y(transmit_i),element_pos_y(mis_index(straight_ray_i)),N);
                    ray_grid_pos_x(1:length(pos_x),mis_index(straight_ray_i),transmit_i) = pos_x';
                    ray_grid_pos_y(1:length(pos_y),mis_index(straight_ray_i),transmit_i) = pos_y';

                    display(['straight ray because of curve failed.  transmiter ', num2str(transmit_i), ', and receiver ', num2str(mis_index(straight_ray_i))]);
                end
            end
        end
    end
      
    % get Coefficient matrix by curve ray tracing result.
    [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    ray tracing for single transmiter with single shoot angle. 
%
%   Input:
%       transmit_pos        transmitter positon (x, y) in meter. coordinate system center locates in the center of recsontruction image. 
%       shoot_angle         shoot angle in Radians
%       ring_radius         transducer ring radius in meter.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
%       refractive_index    n =1540./recon_sound_speed_map; 
%       ray_step_size       grid step size for ray tracing(unit is meter).
% 
%   Output:
%       ray_grid_pos_x:     ray step position x
%       ray_grid_pos_y:     ray step position y
%       tof:                time cost of this ray.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ray_grid_pos_x, ray_grid_pos_y, tof] = single_curve_ray_tracing(transmit_pos, shoot_angle, ring_radius, image_size_p, image_size_m, refractive_index, ray_step_size, speed_in_water)

    % prepare parameter 
    ds = ray_step_size;                                 % ray tracing step size.
    image_resolution = image_size_m/image_size_p;       % resolution of refractive indext map
    current_ray_pos = transmit_pos;                     % initial current ray pos. which starts from transmit position.
    
    % allocate space for save each steps postion and sound distance.
    x = ones(1,3000).*inf;%ˆäãC³1000¨3000
    y = x;
    n_local = zeros(1,3000);%ˆäãC³1000¨3000
    
    % find each steps of curve ray. the longest ray should shortter than 1000 steps
    for ray_step_i = 1:3000 %ˆäãC³1000¨3000
        
        x(ray_step_i) = current_ray_pos(1);         % record current ray position.
        y(ray_step_i) = current_ray_pos(2);
        
        if ray_step_i > 2                           % differential calculation for all steps exceping first step.
            dx = x(ray_step_i)-x(ray_step_i-1);     
            dy = y(ray_step_i)-y(ray_step_i-1);
        else                                        % for the first step.
            dx = ds*cos(shoot_angle);
            dy = ds*sin(shoot_angle);
        end
        
        % current ray position in 
        ray_pos_x_in_image = (image_size_p-1)/image_size_m * (x(ray_step_i)+image_size_m./2) + 1;
        ray_pos_y_in_image = (image_size_p-1)/image_size_m * (y(ray_step_i)+image_size_m./2) + 1;
        
        % first pixel position for bilinear interpolation. nearest refractive map point to current ray pos. 
        rp_x1 = round(ray_pos_x_in_image);
        rp_y1 = round(ray_pos_y_in_image);
        
        
        % differential of refractive index for current ray grid.
        if rp_x1>1 && rp_x1<image_size_p
            nx = (refractive_index(rp_x1+1,rp_y1)-refractive_index(rp_x1-1,rp_y1)) / 2/image_resolution;
        else
            nx = 0;
        end
        if rp_y1>1 && rp_y1<image_size_p
            ny = (refractive_index(rp_x1,rp_y1+1)-refractive_index(rp_x1,rp_y1-1)) / 2/image_resolution;
        else
            ny = 0;
        end
        
        % bilinear interpolation for refractive index of current ray grid positions        
        detx = ray_pos_x_in_image - rp_x1;
        dety = ray_pos_y_in_image - rp_y1;
        
        if detx>=0
            rp_x2 = round(ray_pos_x_in_image+0.5);
        else
            rp_x2 = round(ray_pos_x_in_image-0.5);
        end
        if dety>=0
            rp_y2 = round(ray_pos_y_in_image+0.5);
        else
            rp_y2 = round(ray_pos_y_in_image-0.5);
        end
        
        if rp_x2>0 && rp_x1>0 && rp_x2<=image_size_p && rp_x1<=image_size_p   &&   rp_y1>0 && rp_y2>0 && rp_y1<=image_size_p && rp_y2<=image_size_p
            lx1 = abs(ray_pos_x_in_image-rp_x1);   lx2 = abs(ray_pos_x_in_image-rp_x2);
            ly1 = abs(ray_pos_y_in_image-rp_y1);   ly2 = abs(ray_pos_y_in_image-rp_y2);
            n_inter = refractive_index(rp_x1,rp_y1)*lx2*ly2 + refractive_index(rp_x2,rp_y1)*lx1*ly2 + refractive_index(rp_x1,rp_y2)*lx2*ly1 + refractive_index(rp_x2,rp_y2)*lx1*ly1;
        end
        
        % calculate the next ray grid position
        dsx = (dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*image_resolution^2);
        dsy = (dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*image_resolution^2);
        
        DS = sqrt(dsx^2+dsy^2);   % normaliz to confirm one step is ds
        dsx = dsx/DS*ds;
        dsy = dsy/DS*ds;
        
        current_ray_pos(1) = x(ray_step_i)+dsx;     % next step position x
        current_ray_pos(2) = y(ray_step_i)+dsy;     % next step position y
        n_local(ray_step_i) = n_inter;              % record refractive index of current step.
        
        %‹«ŠE‚ÌÝ’èiŒvŽZ‚ÌI—¹ðŒj boundary condition for ray tracing terminate. 
        if sqrt((current_ray_pos(1)).^2+(current_ray_pos(2)).^2) > ring_radius
            
            % calculate the intersection point between circle and straight line which is defined by last two ray grid position.
            slope_two_point = (current_ray_pos(2)-y(ray_step_i))./(current_ray_pos(1)-x(ray_step_i));           % slope of the straight decided by two point
            
            % get two intersection point between circle ring and the straight.
            if isinf(slope_two_point)
                [intersect_x,intersect_y]=linecirc(slope_two_point , current_ray_pos(1) , 0 , 0 , ring_radius);
            else
                [intersect_x,intersect_y]=linecirc(slope_two_point, current_ray_pos(2)-slope_two_point*current_ray_pos(1) , 0 , 0 , ring_radius);
            end

            % select the intersection point which is close to the last ray grid position.
            if sqrt((intersect_x(1)-current_ray_pos(1))^2+(intersect_y(1)-current_ray_pos(2))^2) < sqrt((intersect_x(2)-current_ray_pos(1))^2+(intersect_y(2)-current_ray_pos(2))^2)     
                x(ray_step_i+1) = intersect_x(1);
                y(ray_step_i+1) = intersect_y(1);
            else
                x(ray_step_i+1) = intersect_x(2);
                y(ray_step_i+1) = intersect_y(2);
            end
            
            break;
        end
    end
    
    % tof of current ray tracing. note: tof has not add the last step time.
    tof = sum(n_local*ds)/speed_in_water;
    step_num = sum(x<inf);
    ray_grid_pos_x = x(1:step_num);
    ray_grid_pos_y = y(1:step_num);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    coefficient projection matrix calculation for straight ray reconstruction method
%
%   Input:
%       ring_diameter       transducer ring diameter in meter.  
%       element_num         transducer elements number.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
% 
%   Output:
%       a_ray_matrix:       ray path weight matrix (coefficient projection matrix) without hanning window.
%       w_ray_matrix:       same to last output but with hanning window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_ray_matrix, w_ray_matrix] = straight_ray_tracing(ring_diameter, element_num, image_size_p, image_size_m)

    
    % preparing the ring and element parameter
    ch = element_num;                                                       % element number
%    theta_rg = linspace(pi/2,(pi/2+(element_num-1)/ch*2*pi),element_num);   % element angles in the ring.
    theta_rg = linspace(0,(element_num-1)/ch*2*pi,element_num);
    
    r_rg = ring_diameter/2;                                                 % radius of the ring.        
    x_rg = r_rg*cos(theta_rg);                                              % element position x
    y_rg = r_rg*sin(theta_rg);                                              % elmeent position y
    
    ray_grid_pos_x = ones(5000,element_num,element_num).*inf;
    ray_grid_pos_y = ones(5000,element_num,element_num).*inf;
    
    
 %   Li = sqrt((element_pos_matrix_x-element_pos_matrix_x')^2+(element_pos_matrix_y-element_pos_matrix_y')^2);     %Œo˜H’·‚³
    
    % prepare recon image parameters
    L = image_size_m/2;             % half length(meter) of recon image
    cell_num = image_size_p;        % length(pixel) of recon image
    cell_size = 2*L/cell_num;       % resolution of recon image
    ray_grid_size = cell_size./2;
    
   
    for transmit_i = 1:element_num
        for receive_i = 1:element_num
            
            Li =  sqrt((x_rg(transmit_i)-x_rg(receive_i))^2+(y_rg(transmit_i)-y_rg(receive_i))^2);  % ray length
            N = round(Li/ray_grid_size);                                                            % step number.
            
            ray_grid_pos_x(1:N,transmit_i, receive_i) = linspace(x_rg(transmit_i),x_rg(receive_i),N);   % position x of each step ????????????? receive transmit position exchanged?????
            ray_grid_pos_y(1:N,transmit_i, receive_i) = linspace(y_rg(transmit_i),y_rg(receive_i),N);   % position y of each step
        end
    end
    
    %  oefficient projection matrix calculation calculation.
    [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m);
                
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    coefficient projection matrix calculation by using positions of each ray step.
%
%   Input:
%       ray_grid_pos_x      position x of each ray steps
%       ray_grid_pos_y      position y of each ray steps.
%       element_num         transducer elements number.
%       image_size_p        image size in pixels 
%       image_size_m        image size in meters
% 
%   Output:
%       a_ray_matrix:       ray path weight matrix (coefficient projection matrix) without hanning window.
%       w_ray_matrix:       same to last output but with hanning window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_ray_matrix, w_ray_matrix] = projection_matrix_calculation(ray_grid_pos_x, ray_grid_pos_y, element_num, image_size_p, image_size_m)
    
%prepare ray beam parameters              
    zz = zeros(image_size_p);           % template for saving coefficients projection matrix for single ray.
    z_cor = zeros(image_size_p);        % template for saving hamming window coefficients projection matrix for single ray. 
    
    % template for save all ray information.
    pos_record_x = zeros(5000,element_num*element_num);     % record ray x position in image.
    zz_record = zeros(5000,element_num*element_num);        % record coefficients projection matrix value.
    z_cor_record = zeros(5000,element_num*element_num);     % record hamming window coefficients projection matrix value.
    
    ray_grid_num = reshape(sum(ray_grid_pos_x<inf),element_num, element_num);
 
    
    % transmit loop
    for transmit_i = 1:element_num   % transmit
        
        % receiver loop
        for receive_i = 1:element_num   %receive
            
            N = ray_grid_num(transmit_i, receive_i);
            x_ry = ray_grid_pos_x(1:N, transmit_i, receive_i);
            y_ry = ray_grid_pos_y(1:N, transmit_i, receive_i);
            
            
            % ray grid loop 
            for kk = 1:N   % ray grid.
                
                %Ä\¬ƒCƒ[ƒW‚ÍƒŠƒ“ƒO’†S‚ð’†S‚Æ‚·‚é‚P•Ól‚Ì³•û—Ìˆæ
                nx2 = (image_size_p-1)/image_size_m * (x_ry(kk)+image_size_m/2)+1;
                ny2 = (image_size_p-1)/image_size_m * (y_ry(kk)+image_size_m/2)+1;
                nx2_1 = round(nx2);         % first pixel x for bilinear interpolation.
                ny2_1 = round(ny2);         % first pixel y for bilinear interpolation.

                %‘‹ŠÖ” hamming window
                w(kk) = 0.54-0.46*cos(2*pi/N*kk);                   % window function, how to define the value for this window.
                
                % find second pixel for bilinear interpolation
                if nx2>nx2_1
                    nx2_2 = round(nx2+0.5);
                else
                    nx2_2 = round(nx2-0.5);
                end
                if ny2>ny2_1
                    ny2_2 = round(ny2+0.5);
                else
                    ny2_2 = round(ny2-0.5);
                end
                
                % conform both pixels are on the recon image.
                if nx2_1>0 && nx2_2>0 && nx2_1<=image_size_p && nx2_2<=image_size_p && ny2_1>0 && ny2_2>0 && ny2_1<=image_size_p && ny2_2<=image_size_p
                    
                    % bilinear interpolation weight calculation.
                    lx1 = abs(nx2-nx2_1);lx2 = abs(nx2-nx2_2);
                    ly1 = abs(ny2-ny2_1);ly2 = abs(ny2-ny2_2);
                    
                    
%                     % note the directions, x is 2nd D direction y is first.
%                     % projection matrix calculation 
%                     zz(ny2_1,nx2_1) = zz(ny2_1,nx2_1)+lx2*ly2;
%                     zz(ny2_1,nx2_2) = zz(ny2_1,nx2_2)+lx1*ly2;
%                     zz(ny2_2,nx2_1) = zz(ny2_2,nx2_1)+lx2*ly1;
%                     zz(ny2_2,nx2_2) = zz(ny2_2,nx2_2)+lx1*ly1;
%                     
%                     % hamming window prjection matrix calculation.
%                     z_cor(ny2_1,nx2_1) = z_cor(ny2_1,nx2_1)+lx2*ly2*w(kk);                    
%                     z_cor(ny2_1,nx2_2) = z_cor(ny2_1,nx2_2)+lx2*ly1*w(kk);
%                     z_cor(ny2_2,nx2_1) = z_cor(ny2_2,nx2_1)+lx1*ly2*w(kk);
%                     z_cor(ny2_2,nx2_2) = z_cor(ny2_2,nx2_2)+lx1*ly1*w(kk);


                    % projection matrix calculation 
                    zz(nx2_1,ny2_1) = zz(nx2_1,ny2_1)+lx2*ly2;
                    zz(nx2_1,ny2_2) = zz(nx2_1,ny2_2)+lx1*ly2;
                    zz(nx2_2,ny2_1) = zz(nx2_2,ny2_1)+lx2*ly1;
                    zz(nx2_2,ny2_2) = zz(nx2_2,ny2_2)+lx1*ly1;
                    
                    % hamming window prjection matrix calculation.
                    z_cor(nx2_1,ny2_1) = z_cor(nx2_1,ny2_1)+lx2*ly2*w(kk);                    
                    z_cor(nx2_1,ny2_2) = z_cor(nx2_1,ny2_2)+lx2*ly1*w(kk);
                    z_cor(nx2_2,ny2_1) = z_cor(nx2_2,ny2_1)+lx1*ly2*w(kk);
                    z_cor(nx2_2,ny2_2) = z_cor(nx2_2,ny2_2)+lx1*ly1*w(kk);

                end
            end
            
            % confirm there is projection ray.
            if sum(zz(:)) >0
                
                pjt_line_num = receive_i + (transmit_i-1)*element_num;  % current projection ray order num.
                
                pos = find(zz>0);               % available pixel position in projection matrix
                size_pos = max(size(pos));      % available pixel size in projection matrix
                pos_record_x(1:size_pos, pjt_line_num) = pos;           % record available pixel position x in projection matrix
                zz_record(1:size_pos,pjt_line_num) = zz(pos);           % record available pixel value in projection matrix
                z_cor_record(1:size_pos,pjt_line_num) = z_cor(pos);     % record available pixel value in hamming window projection matrix
            
                zz(:) = 0;      % seting one ray template to 0
                z_cor(:) = 0;   % seting one ray template to 0;
            end
        end        
    end
    
    % geting whole projection matrix and hamming window projection matrix
    available_mark = find(pos_record_x>0);          % available data position
    xx = pos_record_x(available_mark);              % available data position in its single ray map
    pos_record_y = repmat(1:(element_num^2),[5000,1]);       % data position for ray order number.
    yy = pos_record_y(available_mark);              % available data position for ray order number
    zz_value = zz_record(available_mark);           % coefficients projection matrix data
    z_cor_value = z_cor_record(available_mark);     % hamming window coefficients projection matrix data.
    
    image_resolution = image_size_m/image_size_p;
    ray_grid_resolution = sqrt((ray_grid_pos_x(2,element_num/2,1)-ray_grid_pos_x(1,element_num/2,1))^2  +  (ray_grid_pos_y(2,element_num/2,1)-ray_grid_pos_y(1,element_num/2,1))^2);
    ratio_resolution = ray_grid_resolution./image_resolution;
    
    a_ray_matrix = sparse(xx,yy, zz_value.*ratio_resolution, image_size_p^2, element_num^2);     % coefficients projection matrix sparse matrix
    w_ray_matrix = sparse(xx,yy, z_cor_value.*ratio_resolution, image_size_p^2, element_num^2);  % hamming window coefficients projection matrix sparse matrix.

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    SART reconstruction using coefficient projection matrix calculation
%
%   Input:
%       pjt_data:                       travel time difference map. 
%                                       1st D. is receiver, 
%                                       2nd D. is emitter, 
%       a                               ray path weight matrix (coefficient projection matrix) without hanning window. same to a_ray_matrix in last function.
%       w                               same to last input but with hanning window. same to w_ray_matrix in last function.
%       initial_image                   the departure points of iteration reconstruction (SART).
%       image_size_p                    image size in pixels 
%       image_size_m                    image size in meters
%       iteration_num                   iteration number for calculation.
%       relaxation_factor               relaxation factor for each iteration.
%       available_angle_range           available angle limited during SART
%       speed_in_water                  sound speed in water.
%       max_diff_pjt_emit_recei_pair    pjt pair difference for remove failed pjt data.
%       min_recon_speed                 
%
%   Output:
%       recon_image:            reconstructed ultrasound sound speed image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recon_image = get_recon_image_SART(pjt_data, a, w, initial_image, image_size_p, image_size_m,iteration_num, relaxation_factor, available_angle_range, speed_in_water, max_diff_pjt_emit_recei_pair, min_recon_speed, max_recon_speed)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% mark available data in pjt_data.
    
    % 1.  preparing parameters for angle limited SART
    element_num = size(pjt_data,1);                                             % transducer element number
    image_resolution = image_size_m/image_size_p;                               % image resolution.
    available_receiver_num = round(element_num * available_angle_range/(2*pi)); % available receiver number
    abandoned_receiver_num = element_num - available_receiver_num;              % abandoned receiver numnber which is outof available angle limited. 
    
    %make ang_limited_map
    ang_limited_map = zeros(element_num, element_num);
    temp = cat(1, zeros(round(abandoned_receiver_num/2),1), ones(available_receiver_num,1),zeros((abandoned_receiver_num-round(abandoned_receiver_num/2)),1));
    for i = 1:element_num   % loop for transmiter
        ang_limited_map(:,i) = circshift(temp, [(i-1),0]);
    end
    
    
    % 2. mark failed points in travle time difference map (TTDM) using thereciprocal transmitter?receiver pair.
    fail_limited_map1 = abs(pjt_data - pjt_data') < max_diff_pjt_emit_recei_pair;
    % fail_limited_map1 = ones(size(pjt_data));
    
    
    % 3 limit the range of reconstructed region.
    ring_diameter = 0.2; % ring transducer diameter is 0.1 m
    % t max and t min
    t_max = ring_diameter/min_recon_speed - ring_diameter/speed_in_water;
    t_min = ring_diameter/max_recon_speed - ring_diameter/speed_in_water;
    fail_limited_map2 = pjt_data > t_min & pjt_data< t_max;
    %fail_limited_map2 = ones(size(pjt_data));
    
    
    % 4. combine the unavailable points in TTDM
    position_array = repmat((1:element_num)',[1,element_num]) + (element_num * repmat((0:element_num-1), [element_num, 1]));
    available_pjt_record = position_array(find(ang_limited_map.*fail_limited_map1.*fail_limited_map2));
    

    % get a, w, pjt whith available angle limited.
    a_limited = a(:,available_pjt_record);
    w_limited = w(:,available_pjt_record);
    pjt_limited = pjt_data(available_pjt_record);
    
    % unit exchange for reconstruction. 
    image_g = (image_resolution - (image_resolution/speed_in_water).*initial_image)./initial_image;
    
    % SART iteration
    for iteration_i = 1: iteration_num
        % SART
        image_g(:) = image_g(:) + ((w_limited*((pjt_limited(:) - ((a_limited')*image_g(:)))./sum(a_limited)')) ./ (sum(a_limited,2)+1)).*relaxation_factor;
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing for output and show result for each time iteration. testing
%         recon_image = image_resolution./(image_g+image_resolution/1491);            % unit exchange. 
%         figure; imagesc(recon_image,[1400 1600]);colorbar;                          % show 
%         save(['Method_S_Iteration_',num2str(iteration_i),'.mat'],'recon_image');    % save
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% testing. 
    end
    
    recon_image = image_resolution./(image_g+image_resolution/speed_in_water);    % unit exchange. change unit to sound speed.
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    this function come from matlab "mapping Toolbox" and it is for getting Intersections of circles and lines in Cartesian plane.
%   
%   [xout,yout] = linecirc(slope,intercpt,centerx,centery,radius) finds the points of intersection given a circle defined 
%   by a center and radius in x-y coordinates, and a line defined by slope and y-intercept, or a slope of "inf" and an x-intercept. 
%   Two points are returned. When the objects do not intersect, NaNs are returned.
% 
% When the line is tangent to the circle, two identical points are returned. All inputs must be scalars.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=linecirc(slope,intercpt,centerx,centery,radius)
%LINECIRC  Intersections of circles and lines in Cartesian plane
%
%  [xout,yout] = LINECIRC(slope,intercpt,centerx,centery,radius) finds
%  the points of intersection given a circle defined by a center and
%  radius in x-y coordinates, and a line defined by slope and
%  y-intercept, or a slope of "inf" and an x-intercept.  Two points
%  are returned.  When the objects do not intersect, NaNs are returned.
%  When the line is tangent to the circle, two identical points are
%  returned. All inputs must be scalars
%
%  See also CIRCCIRC.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/11/26 20:35:30 $
% Written by:  E. Brown, E. Byrns

    assert(isscalar(slope) && isscalar(intercpt) && ...
        isscalar(centerx) && isscalar(centery) && isscalar(radius),...
        ['map:' mfilename ':mapError'], 'Inputs must be scalars')

    assert(isreal([slope intercpt centerx centery radius]), ...
        ['map:' mfilename ':mapError'], 'inputs must be real')

    assert(radius > 0, ...
        ['map:' mfilename ':mapError'], 'radius must be positive')

    % find the cases of infinite slope and handle them separately

    if ~isinf(slope)
        % From the law of cosines

        a=1+slope.^2;
        b=2*(slope.*(intercpt-centery)-centerx);
        c=centery.^2+centerx.^2+intercpt.^2-2*centery.*intercpt-radius.^2;

        x=roots([a,b,c])';

        %  Make NaN's if they don't intersect.

        if ~isreal(x)
            x=[NaN NaN]; y=[NaN NaN];
        else
            y=[intercpt intercpt]+[slope slope].*x';
        end

        % vertical slope case
    elseif abs(centerx-intercpt)>radius  % They don't intercept
        x=[NaN;NaN]; y=[NaN;NaN];
    else
        x=[intercpt intercpt];
        step=sqrt(radius^2-(intercpt-centerx)^2);
        y=centery+[step,-step];
    end
end 
