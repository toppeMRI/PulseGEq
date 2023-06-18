function wav = sub_trap2shape(grad, raster)
% Convert trapezoid to (positive) shape on specified raster sample boundary (default: 4us)
% Increase rise/fall times by one sample to keep within slew limit (due to 4us sampling)
%
% Example:
%  seq = mr.Sequence();
%  seq.read('myseq.seq');
%  trap = seq.getBlock(2).gz
%  wav = sub_trap2shape(trap, 4e-6);

if nargin < 2
	raster = 4e-6;   % sec
end

gamma = 4.2576e3;       % Hz/G
gamp = abs(grad.amplitude)/100/gamma;                % Gauss/cm

wav = [ linspace(0, 0, floor(grad.delay/raster) - 1) linspace(0, gamp, ceil(grad.riseTime/raster))  ...
		gamp*ones(1, floor(grad.flatTime/raster)) ... 
		linspace(gamp, 0, ceil(grad.fallTime/raster)) ]';

