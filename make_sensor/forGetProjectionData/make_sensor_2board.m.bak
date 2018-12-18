function [ Sensor ] = make_sensor_2board( Sensor )
%MAKE_SENSOR_2BOARD：2平板トランスデューサのセンサ配列作成関数

Sensor.pos = zeros(2, Sensor.num);%センサ位置[m]
Sensor.pos(1,1:Sensor.num/2) = -Sensor.sizeTotal/2:Sensor.sizeTotal/(Sensor.num/2-1):Sensor.sizeTotal/2 ;%素子水平方向距離[m]
Sensor.pos(2,1:Sensor.num/2) = Sensor.sizeTotal/2;
Sensor.pos(1,Sensor.num/2+1:Sensor.num) = Sensor.pos(1,1:Sensor.num/2);
Sensor.pos(2,Sensor.num/2+1:Sensor.num) = -Sensor.sizeTotal/2;

end

