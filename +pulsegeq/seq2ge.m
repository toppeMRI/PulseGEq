function [modules loopEntries] = seq2ge(seqarg, systemGE, varargin)
% function [modules loopEntries] = seq2ge(seqarg, systemGE, varargin)
%
% Convert a Pulseq file (http://pulseq.github.io/) to a set of TOPPE files
% that can be executed on GE MR scanners. 
%
% See https://toppemri.github.io/ for more info on TOPPE.
%
% This script writes the following files to disk
%   seqstamp.txt
%   *.mod:            One .mod file corresponds to one "unique" block (see below)
%   modules.txt       List of .mod files, and flags indicating whether each .mod file 
%                     corresponds to an RF/ADC/(gradients only) module
%   scanloop.txt      Sequence of instructions for the entire scan (waveform amplitudes, ADC instructions, etc)
%
% Inputs:
%   seqarg            Either a Pulseq file name, or an mr.Sequence object.
%   systemGE          struct        Contains scanner hardware and TOPPE-specific specs. See +toppe/systemspecs.m
% Input options:
%   toppeVersion      string/int    Default is 'v5'/'5'/5 (default) 
%   verbose           boolean       Default: false
%   debug             boolean       Display detailed info about progress (default: false)
%   tarFile           string        default: 'toppeScanFiles.tar'
%   blockStop         int           end at this block in the .seq file (for testing)
%   nt                 Only step through the first nt rows in scanloop.txt. Default: all rows.
%
% Usage example:
%   >> seq = mr.Sequence();
%   >> seq.read('2DFLASH.seq');
%   >> sys = toppe.systemspecs('maxSlew',200,'slewUnit','T/m/s','maxGrad',50','gradUnit','mT/m');
%   >> seq2ge('2DFLASH.seq', sys, 'toppeVersion', 'v6', 'verbose', true);

%% parse inputs
% Defaults
arg.toppeVersion = 6;
arg.verbose = false;
arg.debug = false;
arg.pulseqVersion = 'v1.4.0';
arg.tarFile = 'toppeScanFiles.tar';
arg.blockStop = [];
arg.ibstart = 1;    % skip the first (ibstart-1) events (for testing)
arg.nt      = [];

% Substitute specified system values as appropriate (from MIRT toolbox)
arg = toppe.utils.vararg_pair(arg, varargin);

arg.toppeVersion = arg.toppeVersion(end);
if ischar(arg.toppeVersion)
    arg.toppeVersion = str2num(arg.toppeVersion);
end

switch arg.pulseqVersion
    case 'v1.2.1'
        nEvents = 6;   % number of events per block (number of columns in .seq file)
    case {'v1.3.0', 'v1.3.1', 'v1.4.0'}
        nEvents = 7;
    otherwise
        error(sprintf('Pulseq version %s is not supported', arg.pulseqVersion));
end

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

%% Loop through blocks and build 'modules' and 'loopEntries'

% 'modules' struct array
% Find blocks that are unique in terms of waveforms and timing 
% (i.e., waveform amplitudes, RF/ADC phase, etc can differ),
% and fill 'modules' struct array accordingly. 
% Each entry of 'modules' is a struct containing all waveforms 
% belonging to one module (.mod file), and other module info.
% The usage of the word "module" here is consistent with its usage in TOPPE.
% For now, the 'EXT' event ID (last column in event table) marks the beginning
% of a 'block group' -- this information is used by the GE interpreter.

% 'loopEntries' struct array
% Each entry in this array contains information needed to fill out one row of scanloop.txt.

if arg.verbose
    fprintf('Filling modules struct array, and loopEntries array.\n' );
end

% get contents of [BLOCKS] section
blockEvents = cell2mat(seq.blockEvents);
blockEvents = reshape(blockEvents, [nEvents, length(seq.blockEvents)]).'; 

if ~isempty(arg.blockStop)
    blockEvents = blockEvents(1:arg.blockStop, :);
end

% set number of blocks (rows in .seq file) to step through
if isempty(arg.nt)
    nt = size(blockEvents, 1);
else
    nt = arg.nt;
end

% First entry in 'modules' struct array
block = seq.getBlock(arg.ibstart);
modules(1) = pulsegeq.sub_block2module(block, arg.ibstart, systemGE, 1);

% First entry in 'loopEntries' struct array (first block is by definition a module)
nextblock = seq.getBlock(arg.ibstart+1);   % needed to set 'textra' in scanloop.txt
loopEntries(1) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, 'mod', 1);

% data frames (in Pfile) are stored using indeces 'slice', 'echo', and 'view' 
sl = 1;
view = 1;
echo = 0; 
adcCount = 0;

for ib = (arg.ibstart+1):nt
    if ~mod(ib, 500)
        for inb = 1:20
            fprintf('\b');
        end
        fprintf('Block %d/%d', ib, size(blockEvents, 1));
    end

    block = seq.getBlock(ib);

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
        continue;  % Done, move on to next block 
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
    modCandidate = pulsegeq.sub_block2module(block, ib, systemGE, length(modules) + 1);

    % Is there an existing module that can be reused (scaled)?
    % Specifically, does one of the existing modules (elements of modules) have 
    % the same length waveform, the same non-empty rf/gx/gy/gz, 
    % and the same value of 'hasADC', as modCandidate?
    isUnique = 1; 
    for ic = 1:length(modules)
        if (modules(ic).res == modCandidate.res ...
            & isempty(modules(ic).rf) == isempty(modCandidate.rf) ...
            & isempty(modules(ic).gx) == isempty(modCandidate.gx) ...
            & isempty(modules(ic).gy) == isempty(modCandidate.gy) ...
            & isempty(modules(ic).gz) == isempty(modCandidate.gz) ...
            & modules(ic).hasRF  == modCandidate.hasRF ...
            & modules(ic).hasADC == modCandidate.hasADC ...
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
        modules(end+1) = modCandidate;
        loopEntries(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, ...
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', length(modules));
        continue; % done, so move on to next block
    end

    % modCandidate is not unique and has the same waveform length as an existing module, 
    % so now check to see if all waveform shapes in modCandidate match those in modules(ic).

    tol = 1e-3;  % Shape is deemed equal if norm(difference) < tol
    ii = 1;
    isSameShape = [];
    for ax = {'rf','gx','gy','gz'};
        ax = cell2mat(ax);
        ch = block.(ax);
        if ~isempty(ch)
            for iwav = 1:modules(ic).npulses
                eval(sprintf('wav1 = modules(ic).%s(:,iwav);', ax));
                eval(sprintf('wav2 = modCandidate.%s;', ax));
                if strcmp(ax, 'rf')
                    % subtract common phase before comparing
                    ntmp = round(length(wav1)/2);
                    wav1 = wav1.*exp(-1i*angle(wav1(ntmp)));
                    wav2 = wav2.*exp(-1i*angle(wav2(ntmp)));
                end
                isSameShape(ii,iwav) = norm(wav1-wav2, 1) < tol; % | norm(wav1+wav2) < tol;
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
        loopEntries(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, ...
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic, 'wavnum', iWavReuse);
    else
        % Found a new set of shapes, so add this waveform set to modules(ic)
        modules(ic) = pulsegeq.sub_updatemodule(modules(ic), block, ib, systemGE);
        loopEntries(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, systemGE, ... 
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic, 'wavnum', modules(ic).npulses);
    end

end

if arg.verbose
    fprintf(' done\n');
else
    fprintf('\n');
end

keyboard

%%
%% Write the TOPPE files
%%

%% First, write each module to a .mod file
if arg.verbose
    fprintf(1, 'Writing .mod files and modules.txt...\n');
end

% write modules.txt header
fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n', length(modules));
fprintf(fid,'fname dur(us) hasRF hasADC trigpos hasGrad\n');

% loop through modules
for ic = 1:length(modules)

    hasADC = modules(ic).hasADC;
    hasRF  = modules(ic).hasRF;

    if hasRF & hasADC
        error('Cannot transmit RF and acquire data in same block. Redesign the .seq file.');
    end

    % waveforms in modules are normalized, so now we need to scale to physical units
    rf = [];
    gx = [];
    gy = [];
    gz = [];
    if hasRF
        for ii = 1:modules(ic).npulses
            rfmax= 0;
            for ib = 1:length(loopEntries)
                if loopEntries(ib).mod == ic & loopEntries(ib).wavnum == ii
                    rfmax = max(loopEntries(ib).rfamp, rfmax);
                end
            end
            rf(:,ii) = rfmax * modules(ic).rf(:,ii);
            RFmax(ic,ii) = rfmax;
            %modules(ic).rf(:,ii) = rfmax * modules(ic).rf(:,ii);
        end
    end

    for ii = 1:modules(ic).npulses
        gxmax = 0;
        gymax = 0;
        gzmax = 0;
        for ib = 1:length(loopEntries)
            if loopEntries(ib).mod == ic & loopEntries(ib).wavnum == ii
                gxmax = max(abs(loopEntries(ib).gxamp), gxmax);
                gymax = max(abs(loopEntries(ib).gyamp), gymax);
                gzmax = max(abs(loopEntries(ib).gzamp), gzmax);
            end
        end
        for ax = {'gx','gy','gz'}
            ax = cell2mat(ax);
            eval(sprintf('wav = modules(ic).%s;', ax));
            if ~isempty(wav)
                eval(sprintf('%s(:,ii) = %smax * modules(ic).%s(:,ii);', ax, ax, ax));
                %eval(sprintf('modules(ic).%s(:,ii) = %smax * modules(ic).%s(:,ii);', ax, ax, ax));
            end
        end
        GXmax(ic,ii) = gxmax;
        GYmax(ic,ii) = gymax;
        GZmax(ic,ii) = gzmax;
        %[ic ii rfmax gxmax gymax gzmax]
    end

    %fprintf('module %d: rfmax: %.3f, gxmax:%.2f, gymax:%.2f, gzmax:%.2f\n', ic, rfmax, gxmax, gymax, gzmax);

    if arg.verbose
        fprintf('\tCreating .mod file number %d... ', ic);
    end

    try
        warning('off');    % don't show message about padding waveforms
        nChop(1) = modules(ic).npre;
        nChop(2) = modules(ic).res - modules(ic).npre - modules(ic).rfres;
        %nChop(1) = nChop(1) + mod(nChop(1), 2);  % make even
        %nChop(2) = nChop(2) - mod(nChop(2), 2);
        toppe.writemod(systemGE, 'ofname', modules(ic).fname, ...
            'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, ...
            'nChop', nChop);
        warning('on');
    catch ME
        error(sprintf('Error in writemod:\n%s', ME.message));
    end
        
    if arg.verbose
        fprintf('success\n');
    end

    % update entry in modules.txt
    fprintf(fid,'%s\t%d\t%d\t%d\t-1\n', modules(ic).fname, round(modules(ic).duration), hasRF, hasADC);    

end
fclose(fid);

if arg.verbose
    fprintf('done. Created %d .mod files. To view, type: \t toppe.plotmod(''all'')\n', ic);
end


%% Write scanloop.txt, which specifies the scan sequence (along with modules.txt and the .mod files).

% load .mod files
mods = toppe.tryread(@toppe.readmodulelistfile, 'modules.txt');

% initialize scanloop.txt file
toppe.write2loop('setup', systemGE, 'version', arg.toppeVersion); 

% write to scanloop.txt, one row at a time
for ib = 1:length(loopEntries)

    if isempty(loopEntries(ib).mod)
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

    iMod = loopEntries(ib).mod;
    iWav = loopEntries(ib).wavnum;

    % RF scaling
    if modules(iMod).hasRF
        RFamplitude = loopEntries(ib).rfamp/RFmax(iMod,iWav);
    end

    % gradient scaling
    if GXmax(iMod,iWav) > 0
        Gamplitude(1) = loopEntries(ib).gxamp/GXmax(iMod,iWav);
    end
    if GYmax(iMod,iWav) > 0
        Gamplitude(2) = loopEntries(ib).gyamp/GYmax(iMod,iWav);
    end
    if GZmax(iMod,iWav) > 0
        Gamplitude(3) = loopEntries(ib).gzamp/GZmax(iMod,iWav);
    end

    RFphase  = loopEntries(ib).rfphs;
    DAQphase = loopEntries(ib).recphs;
    RFspoil  = false;
    RFoffset = loopEntries(ib).rffreq;    % Hz
    slice    = loopEntries(ib).slice;
    echo     = loopEntries(ib).echo + 1;  % write2loop starts indexing at 1
    view     = loopEntries(ib).view;
    Dabmodes = {'off','on'};
    dabmode  = Dabmodes{loopEntries(ib).dabmode+1};
    textra   = loopEntries(ib).textra*1e3;    % msec
    core     = loopEntries(ib).blockGroupID;

    if textra < 0
        textraWarning = true;
        textra = 0;
    else
        textraWarning = false;
    end

    %toppe.write2loop(sprintf('module%d.mod',iMod), ...
    toppe.write2loop(modules(iMod).fname, systemGE, ...
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
        'textra',      textra, ...
        'core', core);

end

if textraWarning
    fprintf(['\nWarning: requested textra < 0, which means that .seq sequence timing ', ...
        'is too tight to be directly converted to TOPPE --', ...
        ' ''textra'' set to zero in one or more scanloop.txt entries.\n']);
end

% close file
toppe.write2loop('finish', systemGE);

if arg.toppeVersion > 5
    % Write cores.txt, which defines the block groups
    blockGroups = [];
    for ie=1:length(loopEntries)
        bgID = loopEntries(ie).blockGroupID;
        modID = loopEntries(ie).mod;
        %[ie bgID modID]
        if bgID > 0
            % start of group (will simply overwrite if already existing)
            blockGroups{bgID} = modID;
            bgIDcurrent = bgID;
        else
            blockGroups{bgIDcurrent} = [blockGroups{bgIDcurrent} modID];
        end
    end
    toppe.writecoresfile(blockGroups);
end

% Write TOPPE .entry file.
% This can be edited by hand as needed after copying to scanner.
for ic = 1:length(modules)
    if modules(ic).hasRF
        b1ScalingFile = modules(ic).fname;
    end
    if modules(ic).hasADC
        readoutFile = modules(ic).fname;
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
if arg.toppeVersion > 5
    system(sprintf('tar rf %s %s', arg.tarFile, 'cores.txt'));
end
for ic = 1:length(modules)
    system(sprintf('tar rf %s %s', arg.tarFile, modules(ic).fname));
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

% clean up
system('rm modules.txt scanloop.txt seqstamp.txt toppeN.entry');
for ic = 1:length(modules)
    system(sprintf('rm %s', modules(ic).fname));
end

if ~arg.verbose
    fprintf('\n');
end

fprintf('Remember to rename one of the .mod files to ''tipdown.mod'', and another to ''readout.mod''\n');

return

%% End of main script

