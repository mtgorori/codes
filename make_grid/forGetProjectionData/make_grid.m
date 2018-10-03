function [ Grid ] = make_grid( Sensor )
%MAKE_GRID
%2018/04/02���_�F2018_03_23-_layered_medium����layered_medium_SART_ROI_10mm.m�̂��߂ɍ���Ă���D
% Sensor.num = 256;
% Sensor.sizeTotal = 50.e-3;
cell_size = Sensor.sizeTotal / (Sensor.num-1);
x_grid = -Sensor.sizeTotal/2: cell_size : Sensor.sizeTotal/2;
y_grid = x_grid;
Grid.x = x_grid;
Grid.y = y_grid;
Grid.size = cell_size;

end

