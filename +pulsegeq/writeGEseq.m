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
    gxScale = -1;   % Gx amplitude scaling, int16 
    gyScale = -1; 
    gzScale = -1; 
    recPhsScale = -1;  % rec phase, int16.

    coreID = Dyn(ib,1);

    parentBlockID = Dyn(ib,2);

    if parentBlockID ~= 0  % ignore delay blocks for now. TODO
        parentBlock = ParentBlocks{parentBlockID};

        if parentBlock.maxrfamp > 0
            rfScale = 2*round(0.5*Dyn(ib,3)/parentBlock.maxrfamp*C.MAXIAMP); % hardware units
        end

        rfPhsScale = 2*round(0.5*Dyn(ib,4)/pi*C.MAXIAMP);    % int16: [-32766 32766] = [-pi pi]

        rfFreq = round(Dyn(ib,5));  % int16. Hz

        if parentBlock.maxgxamp > 0
            gxScale = 2*round(0.5*Dyn(ib,6)/parentBlock.maxgxamp*C.MAXIAMP);
        end
        if parentBlock.maxgyamp > 0
            gyScale = 2*round(0.5*Dyn(ib,7)/parentBlock.maxgyamp*C.MAXIAMP);
        end
        if parentBlock.maxgzamp > 0
            gzScale = 2*round(0.5*Dyn(ib,8)/parentBlock.maxgzamp*C.MAXIAMP);
        end

        if Dyn(ib,9) ~= -99
            phs_rad = angle(exp(1i*Dyn(ib,9)));  % wrap to [-pi pi] range
            recPhsScale= 2*round(0.5*phs_rad/pi*C.MAXIAMP); %  int16: [-32766 32766] = [-pi pi]
        else
            recPhsScale = -1;
        end

        d = [coreID parentBlockID rfScale rfPhsScale rfFreq gxScale gyScale gzScale recPhsScale];
        fwrite(fid, d, 'int16');
    end
end

fclose(fid);

return


