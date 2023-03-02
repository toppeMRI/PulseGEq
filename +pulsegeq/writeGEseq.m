function writeGEseq(fname, systemGE, ParentBlocks, Cores, Dyn)
% function writeGEseq(fname, systemGE, ParentBlocks, Cores, Dyn)
%
% ParentBlocks   cell array of Pulseq blocks 
% Cores          block group definitions 
% Dyn            dynamic scan information (see getdynamics.m)

C = pulsegeq.constants;

fid = fopen(fname, 'w', 'ieee-be');

% Parent blocks
fwrite(fid, length(ParentBlocks), 'int16');  % number of parent blocks
for blockID = 1:length(ParentBlocks)
    pulsegeq.writeblock(fid, blockID, ParentBlocks{blockID}, systemGE);
end

% Cores (block groups)
pulsegeq.writecores(fid, Cores);

%% Dynamics (scan loop)
% Convert to hardware units (int16) suitable for execution on GE
% using the PulseGEq interpreter.

% Write number of events/rows in .seq file as two int16 with base 32766.
% We do it this way because interpreter is already set up to read int16.
nt = size(Dyn,1);
fwrite(fid, floor(nt/C.MAXIAMP), 'int16');
fwrite(fid, mod(nt, C.MAXIAMP), 'int16');

% defaults

for ib = 1:nt
    % Dyn(ib, :) = [coreID parentBlockID rfamp rfphs rffreq gxAmp gyAmp gzAmp recphs];
    coreID = Dyn(ib,1);

    parentBlockID = Dyn(ib,2);
    b = ParentBlocks{parentBlockID};

    if parentBlock.maxrfamp > 0
        rfamp4ge = 2*round(0.5*Dyn(ib,3)/parentBlock.maxrfamp*C.MAXIAMP); % hardware units
    else
        rfamp4ge = 0;
    end

    rfphs = Dyn(ib,4);   % radians
    rfphs4ge = 2*round(0.5*Dyn(ib,4)/pi*C.MAXIAMP);    % [pi, pi] = [-32766 32766]

    rffreq4ge = round(Dyn(ib,5));  % Hz

    if pBlock.maxgxamp > 0
        gxAmp4ge = 2*round(0.5*Dyn(ib,6)/parentBlock.maxgxamp*C.MAXIAMP);
    else
        gxAmp4ge = 0;
    end

    recphs4ge = 2*round(0.5*Dyn(ib,6)/pi*C.MAXIAMP);  % [pi, pi] = [-32766 32766]

    d = [coreID parentBlockID rfamp4ge rfphs4ge rffreq4ge recphs4ge gxAmp4ge gyAmp4ge gzAmp4ge];
    fwrite(fid, Dyn(ib,:), 'int16');
end

fclose(fid);

return


function Dyn = getdynamics(block, parentBlock, parentBlockID, coreID)
% Return vector containing waveform amplitudes, RF/ADC phase, etc,
% for a Pulseq block, in hardware units

C = pulsegeq.constants;

% all members must be defined (even if not used)
rfamp = 0;       % int16
gxAmp = 0;       % int16
gyAmp = 0;
gzAmp = 0;
rfphs  = 0;      % int16. [-32766 32766] = [-pi pi]
rffreq = 0;      % int16, Hz
recphs = 0;      % int16. [-32766 32766] = [-pi pi]

if ~isempty(block.rf)
    rfamp = max(abs(block.rf.signal)); % Hz
    if parentBlock.maxrfamp > 0
        rfamp4ge = 2*round(0.5*Dyn(ib,3)/parentBlock.maxrfamp*C.MAXIAMP); % hardware units
    rfphs = block.rf.phaseOffset;      % radians
    rffreq = block.rf.freqOffset;      % Hz
end
if ~isempty(block.gx)
    if parentBlock.maxgxamp > 0
        gxAmp = 2*round(0.5*block.gx.amplitude/parentBlock.maxgxamp*C.MAXIAMP);
    end
end
if ~isempty(block.gy)
    if parentBlock.maxgyamp > 0
        gyAmp = 2*round(0.5*block.gy.amplitude/parentBlock.maxgyamp*C.MAXIAMP);
    end
end
if ~isempty(block.gz)
    if parentBlock.maxgzamp > 0
        gzAmp = 2*round(0.5*block.gz.amplitude/parentBlock.maxgzamp*C.MAXIAMP);
    end
end
if ~isempty(block.adc)
    recphs = 2*round(0.5*block.adc.phaseOffset/pi*C.MAXIAMP);
end

Dyn = [coreID parentBlockID rfamp rfphs rffreq gxAmp gyAmp gzAmp recphs];
    
    
