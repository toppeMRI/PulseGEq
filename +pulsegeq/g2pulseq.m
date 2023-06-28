function gout = g2pulseq(g,geRasterTime,sys)
% convert gradient from Gauss/cm to Hz/m, and interpolate to sys.gradRasterTime
gamma = 4.2576e3;         % Hz/G
g = g * gamma * 100;   % Hz/m
T = numel(g)*geRasterTime;    % pulse duration
tge = geRasterTime/2:geRasterTime:T; %(T-geRasterTime/2);
t = sys.gradRasterTime/2:sys.gradRasterTime:T;
gout = interp1(tge, g, t, 'linear', 'extrap');
return;
