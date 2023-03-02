function Dyn = getdynamics4ge(b, pBlock, parentBlockID, coreID)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in hardware units suitable for the GE Pulseq
% interpreter (int16).

C = pulsegeq.constants;

% define all members (simplifies file format -- does not impact scan)
rfAmp = 0;       % normalized to signed int16
gxAmp = 0;       % normalized to signed int16
gyAmp = 0;       % normalized to signed int16
gzAmp = 0;       % normalized to signed int16
rfphs  = 0;      % normalized to signed int16 (+32766 = +pi)
rffreq = 0;      % Hz, rounded to int16
recphs = 0;      % normalized to signed int16 (+32766 = +pi)

if ~isempty(b.rf)
    rfAmp = max(abs(b.rf.signal)); % Hz
    rfphs = b.rf.phaseOffset;      % radians
    rffreq = b.rf.freqOffset;      % Hz
end
if ~isempty(b.gx)
    if pBlock.maxgxamp > 0
        gxAmp = 2*round(0.5*b.gx.amplitude/pBlock.maxgxamp*C.MAXIAMP);
    end
end
if ~isempty(b.gy)
    if pBlock.maxgyamp > 0
        gyAmp = 2*round(0.5*b.gy.amplitude/pBlock.maxgyamp*C.MAXIAMP);
    end
end
if ~isempty(b.gz)
    if pBlock.maxgzamp > 0
        gzAmp = 2*round(0.5*b.gz.amplitude/pBlock.maxgzamp*C.MAXIAMP);
    end
end
if ~isempty(b.adc)
    recphs = 2*round(0.5*b.adc.phaseOffset/pi*C.MAXIAMP);
end

Dyn = [coreID parentBlockID rfAmp rfphs rffreq recphs gxAmp gyAmp gzAmp];
    
    
