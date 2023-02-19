function d = getdynamicsettings(b)
% Return struct contain waveform amplitudes, RF/ADC phase, etc
% for a Pulseq block

% defaults
d.rfAmplitude = 0;      % Pulseq units (Hz)
d.gxAmplitude = 0;      % Pulseq units (Hz/m)
d.gyAmplitude = 0;      % Pulseq units (Hz/m)
d.gzAmplitude = 0;      % Pulseq units (Hz/m)
d.rfphs       = 0;      % RF phase offset (rad)
d.rffreq      = 0;      % RF frequency offset (Hz)
d.recphs      = 0;      % receive phase (rad)

if ~isempty(b.rf)
    d.rfAmplitude = max(abs(b.rf.signal)); % Hz
    d.rfphs = b.rf.phaseOffset;
    d.rffreq = b.rf.freqOffset;
end
if ~isempty(b.gx)
    d.gxAmplitude = gx.amplitude;
end
if ~isempty(b.gy)
    d.gyAmplitude = gy.amplitude;
end
if ~isempty(b.gz)
    d.gzAmplitude = gz.amplitude;
end
if ~isempty(b.adc)
    d.recphs = b.adc.phaseOffset;
end
    
    
