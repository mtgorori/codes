function [ Medium, Layer ] = make_layered_medium_vertical( thickness, Grid )
%getProjectionDataのため
%make_layered_medium
%概要説明：筋肉，脂肪，筋肉の三層構造．脂肪層の厚さが線形に増加すると，ROI中の平均音速も線形に増加する．
%目的：このデータを基準にして音速分布再構成や音速推定を行って，真値との差を評価する．
%   詳細説明をここに記述
%　input:Layer.thicknessー脂肪層の厚さ[mm]
%　output:Medium.sound_speedー媒質の音速分布[m/s]
%               Medium.densityー媒質の密度分布[kg/m3]
%               Layer.sound_speed_aveーROI中の平均音速[m/s]
%%初期パラメタ
%音速
v_muscle = 1580;%筋肉の音速[m/s]
v_fat = 1450;%脂肪の音速[m/s]
%密度
den_muscle = 1040;%筋肉の密度[kg/m3]
den_fat = 920;%脂肪の密度[kg/m3]
%%処理部
%音速分布生成，密度分布生成
Medium.sound_speed = v_muscle*ones(length(Grid.x),length(Grid.y));
Medium.density = den_muscle*ones(length(Grid.x),length(Grid.y));
Medium.v0 = v_muscle;%背景
Layer.thickness = thickness*1e-3;%[mm]
for i = 1:length(Grid.x)
    if abs(Grid.x(i))<=Layer.thickness
        Medium.sound_speed(:,i) = v_fat;
        Medium.density(:,i) = den_fat;
    end
end
end