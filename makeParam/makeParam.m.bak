function [ param ] = makeParam( dx,dy,Nx,Ny,isRing )
% USCT-Sim‚Ìparam‚ğ¶¬‚·‚éŠÖ”D
param.grid.dx = dx*1e-3;
param.grid.dy = dy*1e-3;
param.grid.Nx = Nx;
param.grid.Ny = Ny;
param.io.save_movie = 1;
if isRing == 1
    prompt = 'Number of the ring arrays?\n';
    param.ringarray.num_points = input(prompt);
    prompt = 'Radius of the ring?\n';
    param.ringarray.radius = input(prompt);
end
param.t_end = 1.5e-04;%realisticscatter‚Ìê‡D
param.sensor.freq = 40000000;
param.source.waveform.freq = 2000000;
param.source.waveform.wavenum = 1;
param.source.waveform.wavenum_offset = 2;
param.source.waveform.magnitude = 50;
param.source.waveform.type = 'sinusoidal';
end

