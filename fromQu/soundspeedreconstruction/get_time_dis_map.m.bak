
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate transmit time table for each pair of pixle and elements with sound speed image.
% writen by Xiaolei Qu in 2015/05/15
%  
% 
% input:    ring_diameter:          a number; ring transducer inner diameter [m]
%           element_num:            a number; tatolly elements number in ring transducer
%           map_grid_num:           a number; grid number for transmit time table (it is 2D for each elements. and x and y direction has same szie which is map_grid_num)
%           sound_speed:            a number or a 2D image; sound speed image or sound speed [m/s] (it could be a number or a image)
%           sound_speed_image_size: a number; this parameter is required if sound_speed is a image which has same size in x and y direction. this is the lenght of image size [m].
%           method_type:            a string; this parameter is required if sound_speed is a image.
%                                   there are five optional methods. 
%                                   1. FMM:     Fast Marching Methods
%                                   2. HAFMM:   High Accuracy Fast Marching Method
%                                   ↓2016/5/23時点で精度が最も高い
%                                   3. MFMM:    Multistencils Fast Marching Method. 1-3 can be found in paper <M. S. Hassouna and A. A. Farag"Multistencils Fast marching methods: A highly accurate solution to the Eikonal equation on Cartesian domains", IEEE Tansactions on PAMI>.
%                                   4. SR:      Sraight ray with sound speed image (ray step Integration).
%                                   5. BR:      Bent ray with sound speed image (ray step integration).
%
% output: time_dis_map:  record transmit time from each pixel to each elements. its size is grid_num*grid_num*element_num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function time_dis_map = get_time_dis_map(varargin)

    % parameter setting
    if nargin == 0
        ring_diameter = 0.2;        % m
        element_num = 2048;
        map_grid_num = 512;
        sound_speed = 1530;         % m/s
        
        time_dis_map = get_time_dis_map_homo(ring_diameter, element_num, map_grid_num, sound_speed);
        
    elseif nargin == 4
        ring_diameter = varargin{1};
        element_num = varargin{2};
        map_grid_num = varargin{3};
        sound_speed = varargin{4};
        
        % method selection
        if numel(sound_speed)==1
            time_dis_map = get_time_dis_map_homo(ring_diameter, element_num, map_grid_num, sound_speed);
        else
            display('ERROR: sound speed image resolution should be inputed for inhomo case');  % if you inputed 4 parameter, the inputed sound_speed should be a number. 
        end
        
    elseif nargin == 8;
        ring_diameter = varargin{1};
        element_num = varargin{2};
        map_grid_num = varargin{3};
        sound_speed_image = varargin{4};    
        sound_speed_image_size_m = varargin{5};
        method_type = varargin{6};
        element_pos_x = varargin{7};
        element_pos_y = varargin{8};
        
        % chercking parameters
        sound_speed_image(sound_speed_image==0) = 1540;
        
        % method selection
        if numel(sound_speed_image)==1
            display('ERROR: sound speed image should be a mtrix for inhomo case');   % if you inputed 6 parameter, the inputed sound_speed should be a distribution image. 
            
        elseif strcmp(method_type,'MFMM') || strcmp(method_type,'HAFMM') || strcmp(method_type,'FMM') % Fast marching method family
            time_dis_map = get_time_dis_map_inhomo(ring_diameter, element_num, map_grid_num, sound_speed_image,sound_speed_image_size_m,method_type,element_pos_x,element_pos_y);
            
        elseif strcmp(method_type,'SR') % straight ray with sound speed image (ray step Integration).
            time_dis_map = get_time_dis_map_inhomo_SR(ring_diameter, element_num, map_grid_num, sound_speed_image,sound_speed_image_size_m,element_pos_x,element_pos_y);
            
        elseif strcmp(method_type,'BR') % bent ray with sound speed image (ray step integration).
            time_dis_map = get_time_dis_map_inhomo_BR(ring_diameter, element_num, map_grid_num, sound_speed_image,sound_speed_image_size_m,element_pos_x,element_pos_y);
        end
        
    end    
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate transmit time from each pixel to each elements.
% with homogenous assumption
% written by Xiaolei Qu in 2015/05/15 
%
% input:  ring_diameter:      ring transducer inner diameter [m]
%         element_num:        tatolly elements number in ring transducer
%         map_grid_num:       the 2D grid has same number in x and y direction.
%         sound_speed:        sound speed in the homogenous [m/s]
% output: time_dis_map:       record transmit time from each pixel to each elements. its size is grid_num*grid_num*element_num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function time_dis_map = get_time_dis_map_homo(ring_diameter, element_num, map_grid_num, sound_speed)

    % cacluate elements coordinate
    element_angle_pos = linspace(0,2*pi*(element_num-1)/element_num,element_num);
    element_pos_x = ring_diameter/2 * cos(element_angle_pos);
    element_pos_y = ring_diameter/2 * sin(element_angle_pos);
    clear element_angle_pos;
    
    % calculate each pixel coordinate
    grid_resolution = (ring_diameter+0.01)/map_grid_num;
    grid_pos = linspace(-(ring_diameter+0.01)/2+grid_resolution/2, (ring_diameter+0.01)/2-grid_resolution/2, map_grid_num);
    [grid_pos_x, grid_pos_y] = meshgrid(grid_pos,grid_pos);
    grid_pos_y = flipud(grid_pos_y);        % transfer meshgrid coordinate system to general cartesian coordinate system.
    clear grid_pos;

    % calculate flight time from each piexle to each element. 
    time_dis_map = sqrt((repmat(grid_pos_x, [1,1,element_num]) - permute(repmat(element_pos_x,[map_grid_num,1,map_grid_num]),[1,3,2])).^2 +...
        (repmat(grid_pos_y, [1,1,element_num]) - permute(repmat(element_pos_y,[map_grid_num,1,map_grid_num]),[1,3,2])).^2)./ sound_speed;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate transmit time from each pixel to each elements using Fast marching method family.
% with sound speed image correction
% written by Xiaolei Qu in 2015/05/15
% 
% input:  ring_diameter:      ring transducer inner diameter [m]
%         element_num:        tatolly elements number in ring transducer
%         map_grid_num:       the 2D grid has same number in x and y direction.
%         sound_speed:        sound speed map obtained by UCT for inhomogenous [m/s]
%         sound_speed_image_size_m:     sond speed map size in meters. 
%         FMM_method:         FMM method type including FMM, HAFMM and MFMM
% output: time_dis_map:       record transmit time from each pixel to each elements. its size is grid_num*grid_num*element_num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function time_dis_map = get_time_dis_map_inhomo(ring_diameter, element_num, map_grid_num, sound_speed_image,sound_speed_image_size_m,FMM_method,element_pos_x,element_pos_y)
    
    % setting FMM method type parameter
    if strcmp(FMM_method,'FMM')
        usecross = false;
        usesecond = false;
    elseif strcmp(FMM_method,'HAFMM')
        usecross = false;
        usesecond = true;
    elseif strcmp(FMM_method,'MFMM')
        usecross = true;
        usesecond = true;
    end
    
    
    % checking size of sound speed image. 
    [input_speed_size_x, input_speed_size_y] = size (sound_speed_image);
    if input_speed_size_x ~= input_speed_size_y
        display('ERROR: the sound speed image should had same size in x and y direction');
        time_dis_map = -1;
        return;
    end
    
    % allocate space for results
    time_dis_map = zeros(map_grid_num,map_grid_num,element_num);
    
    
    % find current sound speed map which corresponds to grid map and come from interplation of inputed sound speed image.
    % input sound speed image pixel coordinate. 
    sound_speed_image_resolution = sound_speed_image_size_m./input_speed_size_x;
    sound_speed_image_grid_pos = linspace(-sound_speed_image_size_m/2+sound_speed_image_resolution/2, sound_speed_image_size_m/2-sound_speed_image_resolution/2, input_speed_size_x);
    [sound_speed_image_grid_pos_x, sound_speed_image_grid_pos_y] = meshgrid(sound_speed_image_grid_pos,sound_speed_image_grid_pos);         % note: current sound speed grid coordinate is in the meshgrid coordinate system.
    
    
    % grid point pixel coordinate
    grid_resolution = sound_speed_image_size_m/map_grid_num;
    grid_pos = linspace(-sound_speed_image_size_m/2, sound_speed_image_size_m/2, map_grid_num);
    [grid_pos_x, grid_pos_y] = meshgrid(grid_pos,grid_pos);
    
    % obtaining the current speed map by interpolation.
    curr_speed_map = interp2(sound_speed_image_grid_pos_x, sound_speed_image_grid_pos_y,sound_speed_image,grid_pos_x,grid_pos_y,'linear');
    curr_speed_map(isnan(curr_speed_map)) = 1540;                       % out side of interpolation region, was set to be 1540 m/s.
    curr_speed_map = curr_speed_map./(grid_pos(2) - grid_pos(1));       % change unit from m/s to pixel/s
    
    

    % get element coordinate
    element_angle_pos = linspace(0,2*pi*(element_num-1)/element_num,element_num);
%     element_pos_x = ring_diameter/2 * cos(element_angle_pos);
%     element_pos_y = ring_diameter/2 * sin(element_angle_pos);←井上修正
%     引数から受け取るように修正
    
    % get element pos in image  (note: the coordinate transformation between image coordinate and general coordinate systems)
    element_pos_x_in_image = interp1(grid_pos, 1:map_grid_num, element_pos_x,'linear'); 
    element_pos_y_in_image = interp1(grid_pos, 1:map_grid_num, element_pos_y,'linear');
    
    
%     % looping for obtain time distance from all pixel to each elements. 
%     for element_i = 1:element_num
%         start_point = ([element_pos_x_in_image(element_i);element_pos_y_in_image(element_i)]);
%         time_dis_map(:,:,element_i) = msfm2d(curr_speed_map, start_point, usesecond, usecross);
%     end

    
%     matlabpool('open', 3);
%     tic;
    for element_i = 1:element_num
%         start_point =
%         round([element_pos_x_in_image(element_i);element_pos_y_in_image(element_i)]);
%       2016/5/24mexファイルのときはmsfm2d、mファイル内の関数を使うときはmsfm2dmat
        time_dis_map(:,:,element_i) = msfm2d(curr_speed_map, [element_pos_x_in_image(element_i);element_pos_y_in_image(element_i)], usesecond, usecross);
    end
%     toc;
%     matlabpool close;
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FMM code start.....
function [T,Y]=msfm2dmat(F, SourcePoints, usesecond, usecross)
    % modified by xiaolei qu 20150423 for sub pixel start position.
    % This function MSFM2D calculates the shortest distance from a list of
    % points to all other pixels in an image, using the  
    % Multistencil Fast Marching Method (MSFM). This method gives more accurate 
    % distances by using second order derivatives and cross neighbours.
    % 
    % T=msfm2d(F, SourcePoints, UseSecond, UseCross)
    %
    % inputs,
    %   F: The speed image. The speed function must always be larger
    %			than zero (min value 1e-8), otherwise some regions will
    %			never be reached because the time will go to infinity. 
    %   SourcePoints : A list of starting points [2 x N] (distance zero)
    %   UseSecond : Boolean Set to true if not only first but also second 
    %                order derivatives are used (default)
    %   UseCross : Boolean Set to true if also cross neighbours 
    %                are used (default)
    % outputs,
    %   T : Image with distance from SourcePoints to all pixels
    %
    % Note:
    %   Compile the c file "mex msfm2d.c"  for cpu-effective registration
    %
    % Literature : M. Sabry Hassouna et Al. Multistencils Fast Marching 
    %   Methods: A Highly Accurate Solution to the Eikonal Equation on
    %   Cartesian Domains
    %
    % Example,
    %   SourcePoint = [51; 51];
    %   SpeedImage = ones([101 101]);
    %   [X Y] = ndgrid(1:101, 1:101);
    %   T1 = sqrt((X-SourcePoint(1)).^2 + (Y-SourcePoint(2)).^2);
    %
    %   % Run fast marching 1th order, 1th order multi stencil 
    %   % and 2th orde and 2th orde multi stencil
    %
    %   tic; T1_FMM1 = msfm2d(SpeedImage, SourcePoint, false, false); toc;
    %   tic; T1_MSFM1 = msfm2d(SpeedImage, SourcePoint, false, true); toc;
    %   tic; T1_FMM2 = msfm2d(SpeedImage, SourcePoint, true, false); toc;
    %   tic; T1_MSFM2 = msfm2d(SpeedImage, SourcePoint, true, true); toc;
    %
    %   % Show results
    %   fprintf('\nResults with T1 (Matlab)\n');
    %   fprintf('Method   L1        L2        Linf\n');
    %   Results = cellfun(@(x)([mean(abs(T1(:)-x(:))) mean((T1(:)-x(:)).^2) max(abs(T1(:)-x(:)))]), {T1_FMM1(:) T1_MSFM1(:) T1_FMM2(:) T1_MSFM2(:)}, 'UniformOutput',false);
    %   fprintf('FMM1:   %9.5f %9.5f %9.5f\n', Results{1}(1), Results{1}(2), Results{1}(3));
    %   fprintf('MSFM1:  %9.5f %9.5f %9.5f\n', Results{2}(1), Results{2}(2), Results{2}(3));
    %   fprintf('FMM2:   %9.5f %9.5f %9.5f\n', Results{3}(1), Results{3}(2), Results{3}(3));
    %   fprintf('MSFM2:  %9.5f %9.5f %9.5f\n', Results{4}(1), Results{4}(2), Results{4}(3));
    %
    % Example multiple starting points,
    %   SourcePoint=rand(2,100)*255+1;
    %   SpeedImage = ones([256 256]);
    %   tic; T1_MSFM2 = msfm2d(SpeedImage, SourcePoint, true, true); toc;
    %   figure, imshow(T1_MSFM2,[]); colormap(hot(256));
    %
    % Function is written by D.Kroon University of Twente (June 2009)

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

%     SourcePoints=int32(floor(SourcePoints));
    
    SourcePoints_start = zeros(2,4*size(SourcePoints,2));
    
    % set all starting points to distance zero and frozen
    for z=1:size(SourcePoints,2)
        % starting point
        x= SourcePoints(1,z); y=SourcePoints(2,z);
        % Set starting point to frozen and distance to zero
        Frozen(ceil(x),ceil(y))=1;      T(ceil(x),ceil(y))=     sqrt((ceil(x)-x).^2 + (ceil(y)-y).^2)/F(ceil(x),ceil(y));
        Frozen(floor(x),floor(y))=1;    T(floor(x),floor(y))=	sqrt((floor(x)-x).^2 + (floor(y)-y).^2)/F(floor(x),floor(y));
        Frozen(ceil(x),floor(y))=1;     T(ceil(x),floor(y))=	sqrt((ceil(x)-x).^2 + (floor(y)-y).^2)/F(ceil(x),floor(y));
        Frozen(floor(x),ceil(y))=1;     T(floor(x),ceil(y))=    sqrt((floor(x)-x).^2 + (ceil(y)-y).^2)/F(floor(x),ceil(y));
        
        SourcePoints_start(1,z+(z-1)*4)     = ceil(x);      SourcePoints_start(2,z+(z-1)*4)     = ceil(y);
        SourcePoints_start(1,z+1+(z-1)*4)   = floor(x);   	SourcePoints_start(2,z+1+(z-1)*4)   = floor(y);
        SourcePoints_start(1,z+2+(z-1)*4)   = ceil(x);   	SourcePoints_start(2,z+2+(z-1)*4)   = floor(y);
        SourcePoints_start(1,z+3+(z-1)*4)   = floor(x);     SourcePoints_start(2,z+3+(z-1)*4)   = ceil(y);
    end

    % Add all neighbours of the starting points to narrow list
    % inoue 素子が位置するピクセルのすぐ隣のピクセルまでの時間のリストをつくる
    for z=1:size(SourcePoints_start,2)
        % starting point
        x=SourcePoints_start(1,z); 
        y=SourcePoints_start(2,z);
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
                    %inoue その点の時間がすでに計算されている場合に小さいほうを選択する（この分岐には現時点の処理では入ってない）
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

    % Loop through all pixels of the image
    for itt=1:numel(F)
        % Get the pixel from narrow list (boundary list) with smallest
        % distance value and set it to current pixel location
        % 現在埋まっている時間リストの中でもっとも小さいものを選択する
        [t,index]=min(neg_list(1,1:neg_pos));
        if(neg_pos==0), break; end
        x=neg_list(2,index); y=neg_list(3,index);
        Frozen(x,y)=1;
        if (itt==1)
            figure;plot(x,y,'bx');xlim([1,512]);ylim([1,512]);
            hold on
        end
        plot(x,y,'bx')
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
    %inoue フローズンされているものを選択する
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
    % then direct neighbours used in solution    %%% mistake, cross position should have different calculation method by xiaolei. 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FMM code over.....




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate transmit time from each pixel to each elements using straight ray.
% with sound speed image correction
% written by Xiaolei Qu in 2015/06/23
% 
% input:  ring_diameter:      ring transducer inner diameter [m]
%         element_num:        tatolly elements number in ring transducer
%         map_grid_num:       the 2D grid has same number in x and y direction.
%         sound_speed:        sound speed map obtained by UCT for inhomogenous [m/s]
%         sound_speed_image_size_m:     sond speed map size in meters. 
% output: time_dis_map:       record transmit time from each pixel to each elements. its size is grid_num*grid_num*element_num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function time_dis_map = get_time_dis_map_inhomo_SR(ring_diameter, element_num, map_grid_num, sound_speed_image,sound_speed_image_size_m,element_pos_x,element_pos_y)
    
    % cacluate elements coordinate
    element_angle_pos = linspace(0,2*pi*(element_num-1)/element_num,element_num);
%     element_pos_x = ring_diameter/2 * cos(element_angle_pos);
%     element_pos_y = ring_diameter/2 * sin(element_angle_pos);
    clear element_angle_pos;
    
    
    % calculate each pixel coordinate 
    grid_resolution = (ring_diameter+0.01)/map_grid_num;
    grid_pos = linspace(-sound_speed_image_size_m/2, sound_speed_image_size_m/2, map_grid_num);
    [grid_pos_x, grid_pos_y] = meshgrid(grid_pos,grid_pos);
    grid_pos_y = flipud(grid_pos_y);    % transfor meshgrid coordinate to cartesian coordinate.
    clear grid_pos;
    
    % calculate step length
    ray_step_length = grid_resolution;
    
    
    
    % get each pixel positions for sound speed image
    input_speed_size_x = size(sound_speed_image,1);
    sound_speed_image_resolution = sound_speed_image_size_m./input_speed_size_x;
    sound_speed_image_pixel_pos = linspace(-sound_speed_image_size_m/2+sound_speed_image_resolution/2, sound_speed_image_size_m/2-sound_speed_image_resolution/2, input_speed_size_x);
    [speed_img_pixel_pos_x, speed_img_pixel_pos_y] =  meshgrid(sound_speed_image_pixel_pos,sound_speed_image_pixel_pos);
    speed_img_pixel_pos_y = flipud(speed_img_pixel_pos_y);    % coordinate transform
    
    
   
    % setting image frame(border or edge) pixel mask
    frame_mask = zeros(map_grid_num);
    frame_mask(1,:) = 1; frame_mask(end,:) = 1; frame_mask(:,1) = 1; frame_mask(:,end) = 1;
    frame_pixel_pos_xy = zeros(map_grid_num*4-4 , 2);
    frame_pixel_pos_xy(:,1) = grid_pos_x(find(frame_mask));
    frame_pixel_pos_xy(:,2) = grid_pos_y(find(frame_mask));
    
    
    % allocate space for saving result time_dis_map
    time_dis_map = zeros(map_grid_num, map_grid_num, element_num);
    parfor emitter_i = 1:1:element_num
%     for emitter_i = 257:257
        tic;
        % allocate space and reset ray step positon and ray step time dis. 
       	ray_step_pos_xy = zeros(round(map_grid_num.*sqrt(2).*2*sum(frame_mask(:))),2) + inf;
        ray_step_time_dis = zeros(round(map_grid_num.*sqrt(2).*2*sum(frame_mask(:))),1) + inf;
        
        
        % get ray step num for each ray and its accumulate value. 
        sample_step_num = round(sqrt((element_pos_x(emitter_i)-frame_pixel_pos_xy(:,1)).^2 + (element_pos_y(emitter_i)-frame_pixel_pos_xy(:,2)).^2)./ ray_step_length);
        sample_step_num_accumulate = cumsum(sample_step_num);
        
        % obtain each step positions. 
        for ray_i = 1:sum(frame_mask(:))
            
            % calculate start and end index in ray_step_pos_xy for current ray 
            if ray_i == 1
                start_index = 1; 
            else
                start_index = sample_step_num_accumulate(ray_i-1) + 1;
            end
            end_index= sample_step_num_accumulate(ray_i);
            
            % get current ray position and save it to ray_step_pos_xy
            ray_step_pos_xy(start_index:end_index,1) = linspace(element_pos_x(emitter_i), frame_pixel_pos_xy(ray_i,1), sample_step_num(ray_i));
            ray_step_pos_xy(start_index:end_index,2) = linspace(element_pos_y(emitter_i), frame_pixel_pos_xy(ray_i,2), sample_step_num(ray_i));

        end
        
        
        % get ray step point sound speed and save it to ray_step_time_dis.
        ray_step_time_dis(1:sample_step_num_accumulate(end)) = interp2(speed_img_pixel_pos_x, speed_img_pixel_pos_y, sound_speed_image, ray_step_pos_xy(1:sample_step_num_accumulate(end),1), ray_step_pos_xy(1:sample_step_num_accumulate(end),2), 'bilinear',1540);
        for ray_i = 1:sum(frame_mask(:))
            % calculate start and end index in ray_step_pos_xy for current ray 
            if ray_i == 1
                start_index = 1; 
            else
                start_index = sample_step_num_accumulate(ray_i-1) + 1;
            end
            end_index= sample_step_num_accumulate(ray_i);
            
            % current step lenght
            current_step_length = sqrt((ray_step_pos_xy(start_index+1,1)-ray_step_pos_xy(start_index,1)).^2 + (ray_step_pos_xy(start_index+1,2)-ray_step_pos_xy(start_index,2)).^2); 
            
            % get ray step point time dist to emitter and save it to ray_step_time_dis. Note: the start position (emitter positon), time should be zeros.
            ray_step_time_dis(start_index:end_index) = cumsum(current_step_length ./ ray_step_time_dis(start_index:end_index)) - current_step_length/ray_step_time_dis(start_index); 
        end
        
        
        % scan coversion interpolation (nearest method)
%         warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        warning('off', 'all');
        time_dis_map(:,:, emitter_i) = griddata(ray_step_pos_xy(1:sample_step_num_accumulate(end),1),ray_step_pos_xy(1:sample_step_num_accumulate(end),2), ray_step_time_dis(1:sample_step_num_accumulate(end)), grid_pos_x, grid_pos_y,'linear');
        warning('on', 'all');
%         warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        toc;
    end
    

end



function time_dis_map = get_time_dis_map_inhomo_BR(ring_diameter, element_num, map_grid_num, sound_speed_image,sound_speed_image_size_m,element_pos_x,element_pos_y)
    
    % cacluate elements coordinate
    element_angle_pos = linspace(0,2*pi*(element_num-1)/element_num,element_num);
%     element_pos_x = ring_diameter/2 * cos(element_angle_pos);
%     element_pos_y = ring_diameter/2 * sin(element_angle_pos);
    clear element_angle_pos;
    
    
    % calculate each pixel coordinate 
    grid_resolution = (ring_diameter+0.01)/map_grid_num;
    grid_pos = linspace(-sound_speed_image_size_m/2, sound_speed_image_size_m/2, map_grid_num);
    [grid_pos_x, grid_pos_y] = meshgrid(grid_pos,grid_pos);
    grid_pos_y = flipud(grid_pos_y);    % transfor meshgrid coordinate to cartesian coordinate.
    clear grid_pos;
    
    % calculate step length
    ray_step_length = grid_resolution;%1540/1.e6/4; % bent ray tracing step should not too small.
    
    
    
    % get each pixel positions for sound speed image
    input_speed_size_x = size(sound_speed_image,1);
    sound_speed_image_resolution = sound_speed_image_size_m./input_speed_size_x;
    sound_speed_image_pixel_pos = linspace(-sound_speed_image_size_m/2+sound_speed_image_resolution/2, sound_speed_image_size_m/2-sound_speed_image_resolution/2, input_speed_size_x);
    [speed_img_pixel_pos_x, speed_img_pixel_pos_y] =  meshgrid(sound_speed_image_pixel_pos,sound_speed_image_pixel_pos);
    speed_img_pixel_pos_y = flipud(speed_img_pixel_pos_y);    % coordinate transform
    
    

    
    ray_tracing_number = element_num*2;
    shoot_angles = linspace(0,2*pi*(ray_tracing_number-1)/ray_tracing_number,ray_tracing_number); 
    
    % allocate space for saving result time_dis_map
    time_dis_map = zeros(map_grid_num, map_grid_num, element_num);
    for emitter_i = 1:element_num
%     for emitter_i = 257:257
        tic;
        emitter_i
        % allocate space and reset ray step positon and ray step time dis. 
       	ray_step_pos_xy = zeros(round(map_grid_num.*sqrt(2).*ray_tracing_number),2) + inf;
        ray_step_time_dis = zeros(round(map_grid_num.*sqrt(2).*ray_tracing_number),1) + inf;
        sample_step_num = zeros(ray_tracing_number, 1);
        
        % get each ray step positions 
        for ray_i = 1:ray_tracing_number
        
            % calculate start and end index in ray_step_pos_xy for current ray 
            if ray_i == 1
                start_index = 1; 
            else
                start_index = sum(sample_step_num(1:ray_i-1)) + 1;
            end
            
            % single_curve_ray_tracing is in matrix coordinate, the input has to be transformed to be matrix coordinate, output has to be transformed from matrix coordinate. 
            [ray_step_pos_x, ray_step_pos_y, sample_step_num(ray_i)] = single_curve_ray_tracing([-element_pos_y(emitter_i),element_pos_x(emitter_i)], shoot_angles(ray_i), input_speed_size_x, sound_speed_image_size_m, sound_speed_image, ray_step_length);
%             if ray_i==1
%                 figure;plot(ray_step_pos_x, ray_step_pos_y)
%                 hold on
%             else
%                 plot(ray_step_pos_x, ray_step_pos_y)
%             end
            
            end_index= sum(sample_step_num(1:ray_i));
            ray_step_pos_xy(start_index:end_index, :) = cat(2, ray_step_pos_y', -ray_step_pos_x');
            
        end
        
        % get ray step point sound speed and save it to ray_step_time_dis.
        ray_step_time_dis(1:sum(sample_step_num(:))) = interp2(speed_img_pixel_pos_x, speed_img_pixel_pos_y, sound_speed_image, ray_step_pos_xy(1:sum(sample_step_num(:)),1), ray_step_pos_xy(1:sum(sample_step_num(:)),2), 'bilinear',1540);
        for ray_i = 1:ray_tracing_number
            % calculate start and end index in ray_step_pos_xy for current ray 
            if ray_i == 1
                start_index = 1; 
            else
                start_index = sum(sample_step_num(1:ray_i-1)) + 1;
            end
            end_index= sum(sample_step_num(1:ray_i));
            
            % current step lenght
            current_step_length = ray_step_length; 
            
            % get ray step point time dist to emitter and save it to ray_step_time_dis. Note: the start position (emitter positon), time should be zeros.
            ray_step_time_dis(start_index:end_index) = cumsum(current_step_length ./ ray_step_time_dis(start_index:end_index)) - current_step_length/ray_step_time_dis(start_index); 
        end
        
        % scan coversion interpolation (nearest method)
        warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        time_dis_map(:,:, emitter_i) = griddata(ray_step_pos_xy(1:sum(sample_step_num(:)),1),ray_step_pos_xy(1:sum(sample_step_num(:)),2), ray_step_time_dis(1:sum(sample_step_num(:))), grid_pos_x, grid_pos_y,'linear');
        warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        toc;
    end

    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description:    ray tracing for single transmiter with single shoot angle. 
%   note: uncompleted.
%
%   Input:
%       transmit_pos        transmitter positon (x, y) in meter. coordinate system center locates in the center of recsontruction image. 
%       shoot_angle         shoot angle in Radians
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
function [ray_grid_pos_x, ray_grid_pos_y, step_num] = single_curve_ray_tracing(transmit_pos, shoot_angle,  image_size_p, image_size_m, sound_speed_image, ray_step_size)

    % prepare parameter 
    ds = ray_step_size;                                 % ray tracing step size.
    image_resolution = image_size_m/image_size_p;       % resolution of refractive indext map
    current_ray_pos = transmit_pos;                     % initial current ray pos. which starts from transmit position.
    speed_in_water = 1540;
    
    % refractive index
    refractive_index = speed_in_water./sound_speed_image;                 % refractive index map
    h = fspecial('average', [3 3]);                 % average filter kernel
    refractive_index = filter2(h,refractive_index);                               % average smoothed refractive index map.
    
    
    % allocate space for save each steps postion and sound distance.
    x = ones(1,1000).*inf;
    y = x;
%     n_local = zeros(1,1000);
    
    % find each steps of curve ray. the longest ray should shortter than 1000 steps
    for ray_step_i = 1:1000 
        
        x(ray_step_i) = current_ray_pos(1);         % record current ray position.
        y(ray_step_i) = current_ray_pos(2);
        if ray_step_i==500
            x_ = x;
        end
            
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
%         n_local(ray_step_i) = n_inter;              % record refractive index of current step.
        
        %境界の設定（計算の終了条件） boundary condition for ray tracing terminate.
        boundary_min_x = -image_size_m/2 + image_resolution/2 + ds;
        boundary_max_x =  image_size_m/2 - image_resolution/2 - ds;
        boundary_min_y = boundary_min_x;
        boundary_max_y = boundary_max_x;
        
        if (current_ray_pos(1)<=boundary_min_x || current_ray_pos(1)>=boundary_max_x) || (current_ray_pos(2)<=boundary_min_y || current_ray_pos(2)>=boundary_max_y)
            x(ray_step_i) = current_ray_pos(1);         % record current ray position.
            y(ray_step_i) = current_ray_pos(2);
            
            break;
        end
        
    end
    
    % tof of current ray tracing. note: tof has not add the last step time.
%     tof = sum(n_local*ds)/speed_in_water;
    step_num = sum(x<inf);
    ray_grid_pos_x = x(1:step_num);
    ray_grid_pos_y = y(1:step_num);
end