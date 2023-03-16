function Dyn = getdynamics(block, parentBlock, parentBlockID, coreID)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in physical (Pulseq) units.
%
% block          Pulseq block obtained with getBlock()
% parentBlock    Parent block that 'block' is a scaled version of
% parentBlockID  int
% coreID         int

C = pulsegeq.constants;

% define all members (simplifies file format -- does not impact scan)
rfamp = 0;       % Hz
gxamp = 0;       % Hz/m
gyamp = 0;
gzamp = 0;
rfphs  = 0;      % radians
rffreq = 0;      % Hz
recphs = -99;    % radians. Value of -99 = no ADC

if ~isempty(block.rf)
    rfamp = max(abs(block.rf.signal));
    rfphs = block.rf.phaseOffset;
    rffreq = block.rf.freqOffset;
end
if ~isempty(block.gx)
%    if parentBlock.maxgxamp > 0
        gxamp = block.gx.amplitude;
%    end
end
if ~isempty(block.gy)
%    if parentBlock.maxgyamp > 0
        gyamp = block.gy.amplitude;
%    end
end
if ~isempty(block.gz)
%    if parentBlock.maxgzamp > 0
        gzamp = block.gz.amplitude;
%    end
end
if ~isempty(block.adc)
    %recphs = 2*round(0.5*block.adc.phaseOffset/pi*C.MAXIAMP);
    recphs = block.adc.phaseOffset;
end

Dyn = [coreID parentBlockID rfamp rfphs rffreq gxamp gyamp gzamp recphs];
    
    
