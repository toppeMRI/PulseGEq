function dyn = getdynamics4ge(b, pBlock, pBlockID, coreID)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in hardware units suitable for the GE Pulseq
% interpreter.

C = pulsegeq.constants;

% define all members (simplifies file format -- does not impact scan)
d.rfAmp = 0;       % normalized to signed int16
d.gxAmp = 0;       % normalized to signed int16
d.gyAmp = 0;       % normalized to signed int16
d.gzAmp = 0;       % normalized to signed int16
d.rfphs  = 0;      % normalized to signed int16 (+32766 = +pi)
d.rffreq = 0;      % Hz, rounded to int16
d.recphs = 0;      % normalized to signed int16 (+32766 = +pi)

if ~isempty(b.rf)
    amp = max(abs(b.rf.signal)); % Hz
    if pBlock.maxrfamp > 0
        d.rfAmp = 2*round(0.5*amp/pBlock.maxrfamp*C.MAXIAMP);
    end
    d.rfphs = 2*round(0.5*b.rf.phaseOffset/pi*C.MAXIAMP);  % [pi, pi] = [-32766 32766]
    d.rffreq = 2*round(0.5*b.rf.freqOffset);  % Hz
end
if ~isempty(b.gx)
    if pBlock.maxgxamp > 0
        d.gxAmp = 2*round(0.5*b.gx.amplitude/pBlock.maxgxamp*C.MAXIAMP);
    end
end
if ~isempty(b.gy)
    if pBlock.maxgyamp > 0
        d.gyAmp = 2*round(0.5*b.gy.amplitude/pBlock.maxgyamp*C.MAXIAMP);
    end
end
if ~isempty(b.gz)
    if pBlock.maxgzamp > 0
        d.gzAmp = 2*round(0.5*b.gz.amplitude/pBlock.maxgzamp*C.MAXIAMP);
    end
end
if ~isempty(b.adc)
    d.recphs = 2*round(0.5*b.adc.phaseOffset/pi*C.MAXIAMP);
end

dyn = [coreID pBlockID d.rfAmp d.rfphs d.rffreq d.recphs d.gxAmp d.gyAmp d.gzAmp];
    
    
