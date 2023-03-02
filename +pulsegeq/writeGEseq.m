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
rfAmp = 0;       % int16
gxAmp = 0;       % int16
gyAmp = 0;
gzAmp = 0;
rfphs  = 0;      % int16. [-32766 32766] = [-pi pi]
rffreq = 0;      % int16, Hz
recphs = 0;      % int16. [-32766 32766] = [-pi pi]

for ib = 1:nt
    % Dyn(ib, :) = [coreID parentBlockID rfAmp rfphs rffreq gxAmp gyAmp gzAmp recphs];
    coreID = Dyn(ib,1);

    parentBlockID = Dyn(ib,2);
    b = ParentBlocks{parentBlockID};

    if parentBlock.maxrfamp > 0
        rfAmp4ge = 2*round(0.5*Dyn(ib,3)/parentBlock.maxrfamp*C.MAXIAMP); % hardware units
    else
        rfAmp4ge = 0;
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

    d = [coreID parentBlockID rfAmp4ge rfphs4ge rffreq4ge recphs4ge gxAmp4ge gyAmp4ge gzAmp4ge];
    fwrite(fid, Dyn(ib,:), 'int16');
end

fclose(fid);
