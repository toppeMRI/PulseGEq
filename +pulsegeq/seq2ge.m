function seq2ge(seqarg, systemGE, varargin)
% function seq2ge(seqarg, systemGE, varargin)
%
% Convert a Pulseq file (http://pulseq.github.io/) to binary file
% that can be executed on GE MR scanners. 

import pulsegeq.*

C = pulsegeq.constants;

%% parse inputs
% Defaults
arg.verbose = false;
arg.debug = false;
arg.nt      = [];

% Substitute specified system values as appropriate (from MIRT toolbox)
arg = pulsegeq.utils.vararg_pair(arg, varargin);


%% Get seq object
if isa(seqarg, 'char')
    seq = mr.Sequence();
    seq.read(seqarg);
else
    if ~isa(seqarg, 'mr.Sequence')
        error('First argument is not an mr.Sequence object');
    end
    seq = seqarg;
end

% get number of blocks (rows in .seq file)
if isempty(arg.nt)
    nEvents = 7;   % Pulseq 1.4.0 
    blockEvents = cell2mat(seq.blockEvents);
    blockEvents = reshape(blockEvents, [nEvents, length(seq.blockEvents)]).'; 
    nt = size(blockEvents, 1);
else
    nt = arg.nt;
end


%% Get parent blocks and write to file
% parent blocks = unique up to a scaling factor, or phase/frequency offsets.
% Contains waveforms with maximum amplitude across blocks.
% First find unique blocks, then determine and set max amplitudes.

% Get unique blocks
ParentBlocks{1} = seq.getBlock(1);

for ib = 1:nt
    if ~mod(ib, 500) | ib == nt
        for inb = 1:20
            fprintf('\b');
        end
        fprintf('Block %d/%d', ib, nt);
    end

    block = seq.getBlock(ib);

    % delay blocks are handled separately (no parent block assigned)
    if (isempty(block.rf) & isempty(block.gx) & isempty(block.gy) & ...
        isempty(block.gz) & isempty(block.adc))
        continue;
    end

    for ipb = 1:length(ParentBlocks)
        IsSame(ipb) = compareblocks(block, ParentBlocks{ipb});
    end
    if sum(IsSame) == 0
        ParentBlocks{ipb+1} = block;  % found a unique block, so add it to list
    end
end

fprintf('\n');

% Determine the parent of each block.
% Pure delay blocks are assigned parent block ID = 0
ParentBlockIDs = zeros(1,nt);   
for ib = 1:nt
    block = seq.getBlock(ib);
    for ipb = 1:length(ParentBlocks)
        if compareblocks(block, ParentBlocks{ipb});
            ParentBlockIDs(ib) = ipb; break;
        end
    end
end

% Determine max amplitude across blocks
for ipb = 1:length(ParentBlocks)
    ParentBlocks{ipb}.maxrfamp = 0;  % this is ok even if rf = []
    ParentBlocks{ipb}.maxgxamp = 0;  
    ParentBlocks{ipb}.maxgyamp = 0;  
    ParentBlocks{ipb}.maxgzamp = 0;  

    for ib = 1:nt
        if ParentBlockIDs(ib) ~= ipb
            continue; 
        end
        block = seq.getBlock(ib);
        if ~isempty(block.rf)
            ParentBlocks{ipb}.maxrfamp = max(ParentBlocks{ipb}.maxrfamp, max(abs(block.rf.signal)));
        end
        if ~isempty(block.gx)
            ParentBlocks{ipb}.maxgxamp = max(ParentBlocks{ipb}.maxgxamp, abs(block.gx.amplitude));
        end
        if ~isempty(block.gy)
            ParentBlocks{ipb}.maxgyamp = max(ParentBlocks{ipb}.maxgyamp, abs(block.gy.amplitude));
        end
        if ~isempty(block.gz)
            ParentBlocks{ipb}.maxgzamp = max(ParentBlocks{ipb}.maxgzamp, abs(block.gz.amplitude));
        end
    end
end

% Set parent block waveform amplitudes to max
for ipb = 1:length(ParentBlocks)
    b = ParentBlocks{ipb};   % shorthand
    if ~isempty(b.rf)
        b.rf.signal = b.rf.signal/max(abs(b.rf.signal))*b.maxrfamp;
    end
    if ~isempty(b.gx)
        b.gx.amplitude = b.maxgxamp;
    end
    if ~isempty(b.gy)
        b.gy.amplitude = b.maxgyamp;
    end
    if ~isempty(b.gz)
        b.gz.amplitude = b.maxgzamp;
    end
end


%% Get cores (block groups) and write to file
currentCoreID = []; 
CoreIDs = zeros(1,nt);  % keep track of which core each block belongs to
for ib = 1:nt
    if ~mod(ib, 500) | ib == nt
        for inb = 1:20
            fprintf('\b');
        end
        fprintf('Block %d/%d', ib, nt);
    end

    b = seq.getBlock(ib);

    if isfield(b, 'label')
        if ~isempty(currentCoreID)  % wrap up the current core
            Cores{currentCoreID} = [currentCoreID length(blocks) blocks];
        end

        % start new core
        currentCoreID = b.label.value;
        blocks = ParentBlockIDs(ib); 
    else
        % add block to currentCoreID
        blocks = [blocks ParentBlockIDs(ib)];
    end

    CoreIDs(ib) = currentCoreID;
end
fprintf('\n');
    

%% Get dynamic scan information
for ib = 1:nt
    b = seq.getBlock(ib);
    parentBlockID = ParentBlockIDs(ib); 
    if parentBlockID ~= 0
        parentBlock = ParentBlocks{parentBlockID}; 
        d(ib,:) = getdynamics(b, parentBlock, parentBlockID, CoreIDs(ib));
    end
end

%% Write sequence to file

pulsegeq.writeGEseq('out.4ge', systemGE, ParentBlocks, Cores, d);


return


%% EOF






%% Old code
if 1
    % get the next block, used to set textra column in scanloop.txt
    if ib < size(blockEvents,1)
        nextblock = seq.getBlock(ib+1);  
    else
        nextblock = [];
    end
%    if isfield(nextblock, 'trig') 
%        nextblock = [];
%    end

    % Empty (pure delay) blocks are accounted for in 'textra' in the previous row in scanloop.txt
    if isempty(block.rf) & isempty(block.adc) ...
        & isempty(block.gx) & isempty(block.gy) & isempty(block.gz) % ...
        %| isfield(block, 'trig')  % ignore trigger (ext) blocks for now. TODO
    end

    % set slice/echo/view indeces (if block is an acquisition block)
    % view = 1, ..., system.maxView
    % sl   = 1, ..., system.maxSlice
    if ~isempty(block.adc)
        view = mod(adcCount, systemGE.maxView) + 1;
        sl   = floor(adcCount/systemGE.maxView) + 1;
        if sl > systemGE.maxSlice;
            error(sprintf('max number of slices ecxeeded (%d)', systemGE.maxSlice));
        end
        echo = floor(adcCount/(systemGE.maxView*systemGE.maxSlice));
        if echo > systemGE.maxEcho
            error(sprintf('max number of echoes ecxeeded (%d)', systemGE.maxEcho));
        end
        %fprintf('ib: %d, view: %d, sl: %d, echo: %d\n', ib, view, sl, echo);

        adcCount = adcCount+1;
    end

    % create a TOPPE module struct from current Pulseq block
    modCandidate = pulsegeq.sub_block2module(block, ib, systemGE, length(moduleArr) + 1);

    % Is there an existing module that can be reused (scaled)?
    % Specifically, does one of the existing modules (elements of moduleArr) have 
    % the same length waveform, the same non-empty rf/gx/gy/gz, 
    % and the same value of 'hasADC', as modCandidate?
    isUnique = 1; 
    for ic = 1:length(moduleArr)
        if (moduleArr(ic).nt == modCandidate.nt ...
            & isempty(moduleArr(ic).rf) == isempty(modCandidate.rf) ...
            & isempty(moduleArr(ic).gx) == isempty(modCandidate.gx) ...
            & isempty(moduleArr(ic).gy) == isempty(modCandidate.gy) ...
            & isempty(moduleArr(ic).gz) == isempty(modCandidate.gz) ...
            & moduleArr(ic).hasRF  == modCandidate.hasRF ...
            & moduleArr(ic).hasADC == modCandidate.hasADC ...
            )
            isUnique = 0;
            break;   % break out of 'for ic' loop. 'ic' now has the value of a module we'll reuse
        end
    end

    if isUnique
        % We found a unique block, so add it as a new module
        if arg.verbose
            fprintf('\tFound new module at block %d\n', ib);
        end
        moduleArr(end+1) = modCandidate;
        loopStructArr(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, ...
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', length(moduleArr));
    end

    % modCandidate is not unique and has the same waveform length as an existing module, 
    % so now check to see if all waveform shapes in modCandidate match those in moduleArr(ic).

    tol = 1e-3;  % Shape is deemed equal if sum(abs(difference)) < tol
    ii = 1;
    isSameShape = [];
    for ax = {'rf','gx','gy','gz'};
        ax = cell2mat(ax);
        ch = block.(ax);
        if ~isempty(ch)
            for iwav = 1:moduleArr(ic).npulses
                eval(sprintf('wav1 = moduleArr(ic).%s(:,iwav);', ax));
                eval(sprintf('wav2 = modCandidate.%s;', ax));
                if strcmp(ax, 'rf')
                    wav1 = abs(wav1);
                    wav2 = abs(wav2);
                end
                isSameShape(ii,iwav) = norm(wav1-wav2,1) < tol;
            end
            ii = ii + 1;
        end
    end

    res = sum(isSameShape,1) == size(isSameShape,1);
    I = find(res==1);
    if ~isempty(I)
        % We found a set of RF/gradient waveforms in modularArr(ic) with the same shapes as those in modCandidate,
        % so we'll reuse that and set 'mod' and 'wavnum' (waveform array column index) accordingly.
        iWavReuse = I(1);
        loopStructArr(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, ...
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic, 'wavnum', iWavReuse);
    else
        % Found a new set of shapes, so add this waveform set to moduleArr(ic)
        moduleArr(ic) = pulsegeq.sub_updatemodule(moduleArr(ic), block, ib, systemGE);
        loopStructArr(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, ... 
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic, 'wavnum', moduleArr(ic).npulses);
    end

end

if arg.verbose
    fprintf(' done\n');
else
    fprintf('\n');
end


%% Optional: As an intermediate check, we can now display the sequence in still/movie mode
%% Example:
if false
    % still frame
    nstart = 1; nstop = 20;
    [rf,gx,gy,gz] = pulsegeq.sub_plotseq(moduleArr,loopStructArr,nstart,nstop);

    % movie
    nBlocksPerTR = 5;
    pulsegeq.sub_playseq(moduleArr, loopStructArr, nBlocksPerTR, nTRskip, tpause);
   pulsegeq.sub_playseq(modArr, loopArr, nBlocksPerTR);
   pulsegeq.sub_playseq(modArr, loopArr, 5, 'gradMode', 'slew', 'tpause', 0.5);
end


%% Hopefully the sequence looks correct (pulsegeq.sub_playseq()), so now we need to write the
%% TOPPE files.


%% First, write each module to a .mod file
if arg.verbose
    fprintf(1, 'Writing .mod files and modules.txt...\n');
end

% write modules.txt header
fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n', length(moduleArr));
fprintf(fid,'wavfile_name    duration (us)     has_RF?     has_ADC?\n');

% loop through moduleArr
for ic = 1:length(moduleArr)

    hasadc = moduleArr(ic).hasADC;
    hasrf  = moduleArr(ic).hasRF;

    if hasrf & hasadc
        error('Cannot transmit RF and acquire data in same block. Redesign the .seq file.');
    end

    % waveforms in moduleArr are normalized, so now we need to scale to physical units
    rf = [];
    gx = [];
    gy = [];
    gz = [];
    if hasrf
        for ii = 1:moduleArr(ic).npulses
            rfmax= 0;
            for ib = 1:length(loopStructArr)
                if loopStructArr(ib).mod == ic & loopStructArr(ib).wavnum == ii
                    rfmax = max(loopStructArr(ib).rfamp, rfmax);
                end
            end
            rf(:,ii) = rfmax * moduleArr(ic).rf(:,ii);
            RFmax(ic,ii) = rfmax;
            %moduleArr(ic).rf(:,ii) = rfmax * moduleArr(ic).rf(:,ii);
        end
    end

    for ii = 1:moduleArr(ic).npulses
        gxmax = 0;
        gymax = 0;
        gzmax = 0;
        for ib = 1:length(loopStructArr)
            if loopStructArr(ib).mod == ic & loopStructArr(ib).wavnum == ii
                gxmax = max(abs(loopStructArr(ib).gxamp), gxmax);
                gymax = max(abs(loopStructArr(ib).gyamp), gymax);
                gzmax = max(abs(loopStructArr(ib).gzamp), gzmax);
            end
        end
        for ax = {'gx','gy','gz'}
            ax = cell2mat(ax);
            eval(sprintf('wav = moduleArr(ic).%s;', ax));
            if ~isempty(wav)
                eval(sprintf('%s(:,ii) = %smax * moduleArr(ic).%s(:,ii);', ax, ax, ax));
                %eval(sprintf('moduleArr(ic).%s(:,ii) = %smax * moduleArr(ic).%s(:,ii);', ax, ax, ax));
            end
        end
        GXmax(ic,ii) = gxmax;
        GYmax(ic,ii) = gymax;
        GZmax(ic,ii) = gzmax;
        %[ic ii rfmax gxmax gymax gzmax]
    end

    %fprintf('module %d: rfmax: %.3f, gxmax:%.2f, gymax:%.2f, gzmax:%.2f\n', ic, rfmax, gxmax, gymax, gzmax);

    if arg.verbose
        fprintf('Creating .mod file number %d...\n', ic);
    end

    % make sure waveforms start and end at zero, and are on a 4-sample boundary (toppe.writemod requires this)
    channels = {'rf','gx','gy','gz'};
    zeropadWarning = false;
    fourSampleBoundaryWarning = false;
    for ii=1:length(channels)
        eval(sprintf('wav = %s;', channels{ii}));
        if ~isempty(wav)
            [nt npulses] = size(wav);
            if any(wav(1,:) ~= 0)
                wav = [zeros(1,npulses); wav];
                zeropadWarning = true;
            end
            if any(wav(end,:) ~= 0)
                wav = [wav; zeros(1,npulses)];
                ZeropadWarning = true;
            end
            [nt npulses] = size(wav);
            if mod(nt,4)
                wav = [wav; zeros(4-mod(nt,4), npulses)];
                fourSampleBoundaryWarning = true;
            end
        end
        eval(sprintf('%s = wav;', channels{ii}));
    end

    if arg.verbose
        if fourSampleBoundaryWarning
            warning(sprintf('One or more waveforms padded with zero at beginning and/or end (module %d).', ic));
        end
        if zeropadWarning
            warning(sprintf('One or more waveforms extended to 16us (4 sample) boundary (module %d).', ic));
        end
    end     

    try
        warning('off');    % don't show message about padding waveforms
        toppe.writemod(systemGE, 'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', moduleArr(ic).ofname); 
        warning('on');
    catch ME
        error(sprintf('Error in writemod:\n%s', ME.message));
    end
        
    if arg.verbose
        fprintf('success\n');
    end

    % update entry in modules.txt
    fprintf(fid,'%s\t%d\t%d\t%d\t-1\n', moduleArr(ic).ofname, 0, hasrf, hasadc);    
end
fclose(fid);

if arg.verbose
    fprintf('done. Created %d .mod files.\n', ic);
    toppe.plotmod('all');
end


%% Write scanloop.txt, which specifies the scan sequence (along with modules.txt and the .mod files).

% load .mod files
mods = toppe.tryread(@toppe.readmodulelistfile, 'modules.txt');

toppe.write2loop('setup', systemGE, 'version', 6); 

for ib = 1:length(loopStructArr)

    if isempty(loopStructArr(ib).mod)
        % skip delay blocks
        continue;
    end

    % defaults
    Gamplitude      = [1 1 1]';
    waveform        = 1;
    textra          = 0;
    RFamplitude     = 1;
    RFphase         = 0;
    DAQphase        = 0;
    RFspoil         = false;
    RFoffset        = 0;
    slice           = 1;
    echo            = 1;
    view            = 1;
    dabmode         = 'on';
    rot             = 0;     % in-plane gradient rotation angle (radians)

    iMod = loopStructArr(ib).mod;
    iWav = loopStructArr(ib).wavnum;

    % RF scaling
    if moduleArr(iMod).hasRF
        RFamplitude = loopStructArr(ib).rfamp/RFmax(iMod,iWav);
    end

    % gradient scaling
    if GXmax(iMod,iWav) > 0
        Gamplitude(1) = loopStructArr(ib).gxamp/GXmax(iMod,iWav);
    end
    if GYmax(iMod,iWav) > 0
        Gamplitude(2) = loopStructArr(ib).gyamp/GYmax(iMod,iWav);
    end
    if GZmax(iMod,iWav) > 0
        Gamplitude(3) = loopStructArr(ib).gzamp/GZmax(iMod,iWav);
    end

    RFphase  = loopStructArr(ib).rfphs;
    DAQphase = loopStructArr(ib).recphs;
    RFspoil  = false;
    RFoffset = loopStructArr(ib).rffreq;    % Hz
    slice    = loopStructArr(ib).slice;
    echo     = loopStructArr(ib).echo + 1;  % write2loop starts indexing at 1
    view     = loopStructArr(ib).view;
    view     = loopStructArr(ib).view;
    Dabmodes = {'off','on'};
    dabmode  = Dabmodes{loopStructArr(ib).dabmode+1};
    textra   = loopStructArr(ib).textra*1e3;    % msec

    if textra < 0
        textraWarning = true;
        textra = 0;
    else
        textraWarning = false;
    end

    %toppe.write2loop(sprintf('module%d.mod',iMod), ...
    toppe.write2loop(moduleArr(iMod).ofname, systemGE, ...
        'Gamplitude',  Gamplitude, ...
        'waveform',    iWav, ...
        'RFamplitude', RFamplitude, ...
        'RFphase',     RFphase, ...
        'DAQphase',    DAQphase, ...
        'RFoffset',    RFoffset, ...
        'slice',       slice, ...
        'echo',        echo, ...
        'view',        view, ...
        'dabmode',     dabmode, ...
        'textra',      textra);

end

if textraWarning
    fprintf(['\nWarning: requested textra < 0, which means that .seq sequence timing ', ...
        'is too tight to be directly converted to TOPPE --', ...
        ' ''textra'' set to zero in one or more scanloop.txt entries.\n']);
end

toppe.write2loop('finish', systemGE);

% Write cores.txt, which defines the block groups
cores = [];
for ie=1:length(loopStructArr)
    bgID = loopStructArr(ie).coreID;
    modID = loopStructArr(ie).mod;
    if ~isempty(bgID)
        % start of group (will simply overwrite if already existing)
        cores{bgID} = modID;
        bgIDcurrent = bgID;
    else
        cores{bgIDcurrent} = [cores{bgIDcurrent} modID];
    end
end
toppe.writecoresfile(cores);

% Write TOPPE .entry file.
% This can be edited by hand as needed after copying to scanner.
for ic = 1:length(moduleArr)
    if moduleArr(ic).hasRF
        b1ScalingFile = moduleArr(ic).ofname;
    end
    if moduleArr(ic).hasADC
        readoutFile = moduleArr(ic).ofname;
    end
end
toppe.writeentryfile('toppeN.entry', ...
    'filePath', '/usr/g/research/pulseq/v6/seq2ge/', ...
    'b1ScalingFile', b1ScalingFile, ...
    'readoutFile', readoutFile);

% Create 'sequence stamp' file for TOPPE.
% This file is listed in line 6 of toppe0.entry
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', systemGE);

% Put TOPPE files in a .tar file (for convenience)
system(sprintf('tar cf %s toppeN.entry seqstamp.txt modules.txt scanloop.txt', arg.tarFile));
system(sprintf('tar rf %s %s', arg.tarFile, 'cores.txt'));
for ic = 1:length(moduleArr)
    system(sprintf('tar rf %s %s', arg.tarFile, moduleArr(ic).ofname));
end

% clean up
system('rm toppeN.entry seqstamp.txt modules.txt scanloop.txt');
system('rm cores.txt');
for ic = 1:length(moduleArr)
    system(sprintf('rm %s', moduleArr(ic).ofname));
end

if arg.verbose
    fprintf(' done\n');
end

return;


% list archive contents
if arg.verbose
    fprintf('\nCreated %s containing the following files:\n', arg.tarFile);
    system(sprintf('tar tf %s', arg.tarFile));
end

if ~arg.verbose
    fprintf('\n');
end

fprintf('Remember to rename one of the .mod files to ''tipdown.mod'', and another to ''readout.mod''\n');

return

%% End of main script

