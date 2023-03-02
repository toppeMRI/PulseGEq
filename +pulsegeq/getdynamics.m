function Dyn = getdynamics4ge(b, parentBlock, parentBlockID, coreID)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block.

C = pulsegeq.constants;

% define all members (simplifies file format -- does not impact scan)
rfAmp = 0;       % Hz
gxAmp = 0;       % Hz/m
gyAmp = 0;
gzAmp = 0;
rfphs  = 0;      % radians
rffreq = 0;      % Hz
recphs = 0;      % radians

if ~isempty(b.rf)
    rfAmp = max(abs(b.rf.signal)); % Hz
    rfphs = b.rf.phaseOffset;      % radians
    rffreq = b.rf.freqOffset;      % Hz
end
if ~isempty(b.gx)
    if parentBlock.maxgxamp > 0
        gxAmp = 2*round(0.5*b.gx.amplitude/parentBlock.maxgxamp*C.MAXIAMP);
    end
end
if ~isempty(b.gy)
    if parentBlock.maxgyamp > 0
        gyAmp = 2*round(0.5*b.gy.amplitude/parentBlock.maxgyamp*C.MAXIAMP);
    end
end
if ~isempty(b.gz)
    if parentBlock.maxgzamp > 0
        gzAmp = 2*round(0.5*b.gz.amplitude/parentBlock.maxgzamp*C.MAXIAMP);
    end
end
if ~isempty(b.adc)
    recphs = 2*round(0.5*b.adc.phaseOffset/pi*C.MAXIAMP);
end

Dyn = [coreID parentBlockID rfAmp rfphs rffreq gxAmp gyAmp gzAmp recphs];
    
    
