function [ projection ] = getProjectionData( Sensor, Grid, Medium )
%GETPROJECTIONDATA：汎用性を高めた投影データ取得関数．
projection = zeros(Sensor.num,Sensor.num);
p = zeros(2,length(Grid.x)+length(Grid.y));

for ii = 1:Sensor.num
    %送信素子の座標設定(x_tr,y_tr)
    x_tr = Sensor.pos(1,ii);
    y_tr = Sensor.pos(2,ii);
    %     pos_tr = [x_tr, y_tr];
    for jj = 1:Sensor.num
        if (ii == jj) || (1<=ii)&&(ii<=Sensor.num/2) && (1<=jj)&&(jj<=Sensor.num/2) || (Sensor.num/2<ii && Sensor.num/2<jj)
            continue
        else
            %受信素子の座標設定(x_re,y_re)
            x_re = Sensor.pos(1,jj);
            y_re = Sensor.pos(2,jj);
            %             pos_re = [x_re, y_re];
            %素子が張る直線の方程式
            x_line = ((x_tr-x_re) / (y_tr-y_re)) * (Grid.y-y_re) + x_re;
            y_line = ((y_tr-y_re) / (x_tr-x_re)) * (Grid.x-x_re) + y_re;
            %各素子位置をグリッドに当てはめる
            [~,x_tr_index] = min(abs(x_tr - Grid.x));
            [~,y_tr_index] = min(abs(y_tr - Grid.y));
            [~,x_re_index] = min(abs(x_re - Grid.x));
            [~,y_re_index] = min(abs(y_re - Grid.y));
            %積分開始位置，終了位置の決定（長方形領域):経路が張る長方形
            x_start = min(x_tr_index,x_re_index);
            x_end = max(x_tr_index,x_re_index);
            y_start = min(y_tr_index,y_re_index);
            y_end = max(y_tr_index,y_re_index);
            %計算領域中の格子との交点を格納(p)
            p(1,1:length(y_start:y_end)) = x_line(y_start:y_end);
            p(2,1:length(y_start:y_end)) = Grid.y(y_start:y_end);
            p(1,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = Grid.x(x_start:x_end);
            p(2,length(y_start:y_end)+1:length(y_start:y_end)+length(x_start:x_end)) = y_line(x_start:x_end);
            p(:,length(y_start:y_end)+length(x_start:x_end)+1:length(p)) = [];
            rm_ind = p(1,:)>Sensor.sizeTotal/2 | p(1,:)<-Sensor.sizeTotal/2 | p(2,:)>Sensor.sizeTotal/2 | p(2,:)<-Sensor.sizeTotal/2;
            p(:,rm_ind) = [];
            p=rmmissing(p,2);%NaNが出たため，対処．
            %初期座標設定(xの値が最小な点）
            if max(p(1,:)) - min(p(1,:)) < Grid.size % argx(min(p))が正しく判定できないことを防ぐため．
                [~,b] = min(p(2,:));
                p_init = p(:,b);
                while ~isempty(p)
                    %線分が属するセルの検出
                    x_ind = find(abs(p_init(1,1) - Grid.x) == 0, 1);
                    y_ind = find(abs(p_init(2,1) - Grid.y) == 0, 1);
                    if isempty(x_ind)
                        [~,x_ind] = min(abs(p_init(1,1) - Grid.x - Grid.size/2));
                    end
                    if isempty(y_ind)
                        [~,y_ind] = min(abs(p_init(2,1) - Grid.y - Grid.size/2));
                    end
                    %線分のもう一方の端点の検出
                    p(:,b) = [];%端点の除去
                    [~,b] = min(p(2,:));%もう一方の端点のインデックス検出
                    p_neighbor = p(:,b);%もう一方の端点検出
                    %寄与音速算出
                    p_length = norm(p_init - p_neighbor);
                    projection(jj,ii) = projection(jj,ii) + p_length*(1/Medium.sound_speed(y_ind,x_ind) - 1/Medium.v0);
                    p_init = p_neighbor;
                end
            else
                
                [~,b] = min(p(1,:));
                p_init = p(:,b);
                while ~isempty(p)
                    %線分が属するセルの検出
                    x_ind = find(abs(p_init(1,1) - Grid.x) == 0, 1);
                    y_ind = find(abs(p_init(2,1) - Grid.y) == 0, 1);
                    if isempty(x_ind)
                        [~,x_ind] = min(abs(p_init(1,1) - Grid.x - Grid.size/2));
                    end
                    if isempty(y_ind)
                        [~,y_ind] = min(abs(p_init(2,1) - Grid.y - Grid.size/2));
                    end
                    %線分のもう一方の端点の検出
                    p(:,b) = [];%端点の除去
                    [~,b] = min(p(1,:));%もう一方の端点のインデックス検出
                    p_neighbor = p(:,b);%もう一方の端点検出
                    %寄与音速算出
                    p_length = norm(p_init - p_neighbor);
                    projection(jj,ii) = projection(jj,ii) + p_length*(1/Medium.sound_speed(y_ind,x_ind) - 1/Medium.v0);
                    p_init = p_neighbor;
                end
            end
        end
    end
end

end

