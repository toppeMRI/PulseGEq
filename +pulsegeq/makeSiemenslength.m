function g = makeSiemenslength(g)
% function g = makeSiemenslength(g)
%
% If length not divisible by 20, pad with zeroes at end
%
% This is done to ensure that GE gradient (raster time 4us)
% has same duration as Siemens gradient (raster time 10us) after
% interpolation (in ge2seq.m)
%
% Input:
%   g    waveform array [nt ...]
%

nb = 20;
if mod(size(g,1),nb)
	ncols = size(g,2);
	g = [g; zeros(nb-mod(length(g),nb),ncols)];
end

% EOF
