function gout = g2pulseq(g,geRasterTime,sys)
% convert gradient from Gauss/cm to Hz/m, and interpolate to sys.gradRasterTime
gamma = 4.2576e3;         % Hz/G
g = g * gamma * 100;   % Hz/m
T = numel(g)*geRasterTime;    % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:sys.gradRasterTime:(T-sys.gradRasterTime);
gout = interp1(tge, g, t, 'linear', 'extrap');
return;
