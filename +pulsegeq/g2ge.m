function gout = g2ge(g,geRasterTime,seq)
% convert gradient from Gauss/cm to Hz/m, and interpolate to seq.gradRasterTime
gamma = 4.2576e3;         % Hz/G
g = g * gamma * 100;   % Hz/m
T = numel(g)*geRasterTime;    % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.gradRasterTime:(T-seq.gradRasterTime);
gout = interp1(tge, g, t, 'linear', 'extrap');
return;
