function gout = g2ge(g,geRasterTime,seq)
% convert gradient from Hz/m to Gauss/cm, and interpolate to geRasterTime
gamma = 4.2576e3;         % Hz/G
g = g / gamma / 100;   % G/cm
T = numel(g)*seq.gradRasterTime;    % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.gradRasterTime:(T-seq.gradRasterTime);
gout = interp1(t, g, tge, 'linear', 'extrap');
return;
