function [ Grid ] = make_grid( Sensor )
%MAKE_GRID
%2018/04/02時点：2018_03_23-_layered_medium内のlayered_medium_SART_ROI_10mm.mのために作ってある．
% Sensor.num = 256;
% Sensor.sizeTotal = 50.e-3;
cell_size = Sensor.sizeTotal / (Sensor.num-1);
x_grid = -Sensor.sizeTotal/2: cell_size : Sensor.sizeTotal/2;
y_grid = x_grid;
Grid.x = x_grid;
Grid.y = y_grid;
Grid.size = cell_size;

end

