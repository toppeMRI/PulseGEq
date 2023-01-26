function rfout = rf2pulseq(rf,geRasterTime,seq)
% function rfout = rf2pulseq(rf,geRasterTime,seq)
% convert rf units from Gauss to Hz, and interpolate to seq.rfRasterTime
%
% rf             [n 1]  RF waveform, Gauss
% geRasterTime   [1]   sec
% seq            Pulseq object. Only seq.rfRasterTime is used.
gamma = 4.2576e3;       % Hz/G
rf = rf*gamma;       % Hz
T = numel(rf)*geRasterTime;   % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.rfRasterTime:(T-seq.rfRasterTime);
rfout = interp1(tge, rf, t, 'linear', 'extrap');
%L = 10; cutoff = 0.9;
%rf = interp(rf,dt/rfRasterTime,L,cutoff);      % upsample from 4us to 1us
return;

