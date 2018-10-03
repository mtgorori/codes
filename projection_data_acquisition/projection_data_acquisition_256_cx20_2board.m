clear
close all
%%%% 初期パラメータ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%音速
v_water = 1540;%水の音速[m/s]
v_fat = 1420;%脂肪の音速[m/s]
%脂肪領域形状
fat_cx = 20.e-3; fat_cy = 0; % 中心位置
fat_radius = 10.e-3;
fd = fat_radius;
%センサ設置
t_size = 100.e-3;
t_num = 256;%トランスデューサ数
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%素子水平方向距離[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
%セル設定
cell_num = 256;
p_size = 50.e-3;
cell_size = p_size*2 / cell_num;
x_grid = -p_size: cell_size : p_size;
y_grid = x_grid;


%%%% 音速分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_dist = v_water*ones(length(x_grid),length(y_grid));
for i = 1:length(x_grid)
    for j = 1:length(y_grid)
        if (x_grid(i)-fat_cx)^2 + (y_grid(j)-fat_cy)^2 <= fd^2
            v_dist(j,i) = v_fat;
        end
    end
end

%%%% 投影データ取得 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projection = zeros(t_num,t_num);
p = zeros(2,length(x_grid)+length(y_grid));

for ii = 1:t_num
    %送信素子の座標設定(x_tr,y_tr)
    x_tr = t_pos(1,ii);
    y_tr = t_pos(2,ii);
    pos_tr = [x_tr, y_tr];
    for jj = 1:t_num
        if (ii == jj) || (1<=ii)&&(ii<=t_num/2) && (1<=jj)&&(jj<=t_num/2) || (t_num/2<ii && t_num/2<jj)
            continue
        else
            %受信素子の座標設定(x_re,y_re)
            x_re = t_pos(1,jj);
            y_re = t_pos(2,jj);
            pos_re = [x_re, y_re];
            %素子が張る直線の方程式
            x_line = ((x_tr-x_re) / (y_tr-y_re)) * (y_grid-y_re) + x_re;
            y_line = ((y_tr-y_re) / (x_tr-x_re)) * (x_grid-x_re) + y_re;
            %各素子位置をグリッドに当てはめる
            [~,x_tr_index] = min(abs(x_tr - x_grid));
            [~,y_tr_index] = min(abs(y_tr - y_grid));
            [~,x_re_index] = min(abs(x_re - x_grid));
            [~,y_re_index] = min(abs(y_re - y_grid));
            %積分開始位置，終了位置の決定（長方形領域):経路が張る長方形
            x_start = min(x_tr_index,x_re_index);
            x_end = max(x_tr_index,x_re_index);
            y_start = min(y_tr_index,y_re_index);
            y_end = max(y_tr_index,y_re_index);
            %計算領域中の格子との交点を格納(p)
            p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
            p(2,1:length(y_start:y_end)) = y_grid(y_start:y_end);
            p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = x_grid(x_start:x_end);
            p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
            p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
            rm_ind = find(p(1,:)>t_size/2 | p(1,:)<-t_size/2 | p(2,:)>t_size/2 | p(2,:)<-t_size/2);
            p(:,rm_ind) = [];
            p=rmmissing(p,2);%NaNが出たため，対処．
            %初期座標設定(xの値が最小な点）
            if max(p(1,:)) - min(p(1,:)) < cell_size % argx(min(p))が正しく判定できないことを防ぐため．
                [~,b] = min(p(2,:));
                p_init = p(:,b);
                while ~isempty(p)
                    %線分が属するセルの検出
                    x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                    y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                    if isempty(x_ind)
                        [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                    end
                    if isempty(y_ind)
                        [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                    end
                    %線分のもう一方の端点の検出
                    p(:,b) = [];%端点の除去
                    [~,b] = min(p(2,:));%もう一方の端点のインデックス検出
                    p_neighbor = p(:,b);%もう一方の端点検出
                    %寄与音速算出
                    p_length = norm(p_init - p_neighbor);
                    projection(jj,ii) = projection(jj,ii) + p_length*(1/v_dist(y_ind,x_ind) - 1/v_water);
                    p_init = p_neighbor;
                end
            else
                
                [~,b] = min(p(1,:));
                p_init = p(:,b);
                while ~isempty(p)
                    %線分が属するセルの検出
                    x_ind = find(abs(p_init(1,1) - x_grid) == 0, 1);
                    y_ind = find(abs(p_init(2,1) - y_grid) == 0, 1);
                    if isempty(x_ind)
                        [~,x_ind] = min(abs(p_init(1,1) - x_grid - cell_size/2));
                    end
                    if isempty(y_ind)
                        [~,y_ind] = min(abs(p_init(2,1) - y_grid - cell_size/2));
                    end
                    %線分のもう一方の端点の検出
                    p(:,b) = [];%端点の除去
                    [~,b] = min(p(1,:));%もう一方の端点のインデックス検出
                    p_neighbor = p(:,b);%もう一方の端点検出
                    %寄与音速算出
                    p_length = norm(p_init - p_neighbor);
                    projection(jj,ii) = projection(jj,ii) + p_length*(1/v_dist(y_ind,x_ind) - 1/v_water);
                    p_init = p_neighbor;
                end
            end
        end
    end
end

%%%% 描画用その1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%CT形状
axis square
hold on
plot(t_pos(1,:),t_pos(2,:),'ro',...
    'MarkerSize',5)
%脂肪領域形状
imagesc(x_grid,y_grid,v_dist)
hold off
%%%% 描画用その2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc(projection);
colorbar;
set(gca,'YDir','normal');
xlabel('Transmitter')
ylabel('Receiver')