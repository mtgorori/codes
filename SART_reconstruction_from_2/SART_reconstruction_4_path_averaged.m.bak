%%%% 初期パラメータ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CT形状
ring_cx = 0; ring_cy = 0;% 中心位置
ring_radius = 50.e-3;% 半径
rr =  ring_radius;
%音速
v_water = 1540;%水の音速[m/s]
v_fat = 1420;%脂肪の音速[m/s]
%脂肪領域形状
fat_cx = 0; fat_cy = 0; % 中心位置
fat_radius = 10.e-3;
fd = fat_radius;
%センサ設置
t_num = 256;%トランスデューサ数
t = linspace(0,2*pi,t_num);%センサ位置角
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,:) = rr*cos(t)+ring_cx;
t_pos(2,:) = rr*sin(t)+ring_cy;
%セル設定
cell_num = 256;
cell_size = ring_radius*2 / cell_num;
x_grid = -ring_radius: cell_size : ring_radius;
y_grid = x_grid;
[X,Y] = meshgrid(x_grid,y_grid);%x,yのメッシュ作成

%%%% 音速分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_dist = v_water*ones(length(x_grid),length(y_grid));
for i = 1:length(x_grid)
    for j = 1:length(y_grid)
        if (x_grid(i)-fat_cx)^2 + (y_grid(j)-fat_cy)^2 <= fd^2
            v_dist(i,j) = v_fat;
        end
    end
end

%%%% 音速再構成 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_sim = zeros(length(x_grid),length(y_grid));
ray_dist = zeros(length(x_grid),length(y_grid));
p = zeros(2,length(x_grid)+length(y_grid));

for ii = 1:t_num
    %送信素子の座標設定(x_tr,y_tr)
    x_tr = t_pos(1,ii);
    y_tr = t_pos(2,ii);
    pos_tr = [x_tr, y_tr];
    for jj = 1:t_num
        if ii == jj
            continue
        else
            v_pre_ave = 0;%平均音速分布算出準備のための変数．
            %受信素子の座標設定(x_re,y_re)
            x_re = t_pos(1,jj);
            y_re = t_pos(2,jj);
            pos_re = [x_re, y_re];
            all_length = norm(pos_tr - pos_re);
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
            rm_ind = find(sqrt(p(1,:).^2 + p(2,:).^2)>rr);
            p(:,rm_ind) = [];
            %初期座標設定(xの値が最小な点）
            [~,b] = min(p(1,:));
            p_init = p(:,b);
            q = zeros(3,length(p));
            q_count = 1;
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
                v_pre_ave = v_pre_ave + v_dist(x_ind,y_ind) * (p_length / all_length);
                ray_dist(x_ind,y_ind) = ray_dist(x_ind,y_ind) +  p_length / all_length;
                p_init = p_neighbor;
                %一時保存
                q(1,q_count) = p_length;
                q(2,q_count) = x_ind;
                q(3,q_count) = y_ind;
                q_count =  q_count + 1;
            end
            for kk = 1:length(q)
                v_sim(q(2,kk),q(3,kk)) = v_pre_ave * q(1,kk) / all_length;
            end
        end
    end
end


v_sim = v_sim ./ ray_dist;

%%%% 描画用その1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%CT形状
plot(t_pos(1,:),t_pos(2,:))
axis square
hold on
plot(t_pos(1,:),t_pos(2,:),'ro',...
    'MarkerSize',1)
%脂肪領域形状
fill(fd*cos(t)+fat_cx,fd*sin(t)+fat_cy,'y')
hold off
%%%% 描画用その2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc(x_grid,y_grid,v_sim);