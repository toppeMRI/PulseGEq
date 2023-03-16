function writeGEseq(fname, systemGE, ParentBlocks, Cores, Dyn)
% function writeGEseq(fname, systemGE, ParentBlocks, Cores, Dyn)
%
% Write sequence to a binary file for the GE Pulseq interpreter.
%
% systemGE       struct containing scanner system info, see systemspecs.m
% ParentBlocks   cell array of Pulseq blocks 
% Cores          block group definitions 
% Dyn            dynamic scan information, in physical units (see getdynamics.m)

C = pulsegeq.constants;

fid = fopen(fname, 'w', 'ieee-be');

%% Parent blocks
fwrite(fid, length(ParentBlocks), 'int16');  % number of parent blocks
for blockID = 1:length(ParentBlocks)
    pulsegeq.writeblock(fid, blockID, ParentBlocks{blockID}, systemGE);
end

%% Cores (block groups)
pulsegeq.writecores(fid, Cores);

%% Dynamics (scan loop)
% Convert to hardware units (int16) suitable for execution on GE
% using the PulseGEq interpreter.

% Write number of events/rows in .seq file as two int16 with base 32766.
% We do it this way because interpreter is already set up to read int16.
nt = size(Dyn,1);
fwrite(fid, floor(nt/C.MAXIAMP), 'int16');
fwrite(fid, mod(nt, C.MAXIAMP), 'int16');

for ib = 1:nt
    % Defaults. A value of -1 tells the interpreter that no pulse exists
    % for that entry.
    rfScale = -1;   % RF amplitude scaling, int16
    gxTrapScale = -1;   % Gx trapezoid amplitude scaling, int16 
    gyTrapScale = -1;
    gzTrapScale = -1;
    gxArbScale = -1;       % Gx arbitrary waveform amplitude scaling, int16 
    gyArbScale = -1; 
    gzArbScale = -1; 

    coreID = Dyn(ib,1);

    parentBlockID = Dyn(ib,2);

    if parentBlockID ~= 0  % ignore delay blocks for now. TODO
        parentBlock = ParentBlocks{parentBlockID};

        if parentBlock.maxrfamp > 0
            rfScale = 2*round(0.5*Dyn(ib,3)/parentBlock.maxrfamp*C.MAXIAMP); % hardware units
        end

        phs_rad = angle(exp(1i*Dyn(ib,4)));  % wrap to [-pi pi] range
        rfPhsScale = 2*round(0.5*phs_rad/pi*C.MAXIAMP);    % int16: [-32766 32766] = [-pi pi]

        rfFreq = round(Dyn(ib,5));  % int16. Hz

        if ~isempty(parentBlock.gx)
            if strcmp(parentBlock.gx.type, 'trap')
                gxTrapScale = 2*round(0.5*Dyn(ib,6)/parentBlock.maxgxamp*C.MAXIAMP);
            else
                gxArbScale = 2*round(0.5*Dyn(ib,6)/parentBlock.maxgxamp*C.MAXIAMP);
            end
        end
        if ~isempty(parentBlock.gy)
            if strcmp(parentBlock.gy.type, 'trap')
                gyTrapScale = 2*round(0.5*Dyn(ib,7)/parentBlock.maxgyamp*C.MAXIAMP);
            else
                gyArbScale = 2*round(0.5*Dyn(ib,7)/parentBlock.maxgyamp*C.MAXIAMP);
            end
        end
        if ~isempty(parentBlock.gz)
            if strcmp(parentBlock.gz.type, 'trap')
                gzTrapScale = 2*round(0.5*Dyn(ib,8)/parentBlock.maxgzamp*C.MAXIAMP);
            else
                gzArbScale = 2*round(0.5*Dyn(ib,8)/parentBlock.maxgzamp*C.MAXIAMP);
            end
        end

        if Dyn(ib,9) ~= -99
            phs_rad = angle(exp(1i*Dyn(ib,9)));  % wrap to [-pi pi] range
            recPhsScale= 2*round(0.5*phs_rad/pi*C.MAXIAMP); %  int16: [-32766 32766] = [-pi pi]
        else
            recPhsScale = -1;
        end

        Dyn_hw = [coreID parentBlockID rfScale rfPhsScale rfFreq gxTrapScale gyTrapScale gzTrapScale gxArbScale gyArbScale gzArbScale recPhsScale];
        fwrite(fid, Dyn_hw, 'int16');
    end
end

fclose(fid);

return


