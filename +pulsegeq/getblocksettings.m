function d = getblocksettings(b)
%
% Return struct contain waveform amplitudes, phase, etc

d.rfAmplitude = 0;      % Pulseq units (Hz)
d.gxAmplitude = 0;      % Pulseq units (Hz/m)
d.gyAmplitude = 0;      % Pulseq units (Hz/m)
d.gzAmplitude = 0;      % Pulseq units (Hz/m)

if ~isempty(b.rf)
    
