function rfout = rf2ge(rf,geRasterTime,seq)
% convert rf units from Hz to Gauss, and increase raster time from 1us to 4us

raster = 4;  % us
gamma = 4.2576e3;    % Hz/G
rf = rf/gamma;       % G
rfout = [0; (rf(2:raster:end)+rf(3:raster:end))/2; 0];

return;

