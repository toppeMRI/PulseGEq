function [modules loop] = seq2ge(seqarg, systemGE, varargin)
% function [modules loop] = seq2ge(seqarg, systemGE, varargin)
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
arg.toppeVersion = 5;
arg.verbose = false;
arg.debug = false;
arg.pulseqVersion = 'v1.4.0';
arg.tarFile = 'toppeScanFiles.tar';
arg.blockStop = [];
arg.nstart = 1;    % skip the first (nstart-1) events (for testing)
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

%% Loop through blocks and build 'modules' and 'loop'

% 'modules' struct array
% Find blocks that are unique in terms of waveforms and timing 
% (i.e., waveform amplitudes, RF/ADC phase, etc can differ),
% and fill 'modules' array accordingly. 
% Each element of 'modules' is a struct containing all waveforms 
% belonging to one module (.mod file), and other module info.
% The usage of the word "module" here is consistent with its usage in TOPPE.

% TODO: add support for groups back in. Need to decide how.
% Was: the 'EXT' event ID (last column in event table) marks the beginning
% of a 'block group' -- this information is used by the GE interpreter.

% 'loop' struct array
% Each entry in this array contains information needed to fill out one row of scanloop.txt.

if arg.verbose
    fprintf('Filling modules struct array, and loop array.\n' );
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
block = seq.getBlock(arg.nstart);
modules(1) = pulsegeq.sub_block2module(block, arg.nstart, systemGE, 1);

% First entry in 'loop' struct array (first block is by definition a module)
nextBlock = seq.getBlock(arg.nstart+1);   % needed to set 'textra' in scanloop.txt
loop(1) = pulsegeq.sub_updateloopstruct([], block, nextBlock, systemGE, 'mod', 1);

% data frames (in Pfile) are stored using indeces 'slice', 'echo', and 'view' 
sl = 1;
view = 1;
echo = 0; 
adcCount = 0;

for n = (arg.nstart+1):nt
    if ~mod(n, 500) | n == nt
        for inb = 1:20
            fprintf('\b');
        end
        fprintf('Block %d/%d', n, size(blockEvents, 1));
    end

    block = seq.getBlock(n);

    % get the next block. If pure delay block, use it to set textra column in scanloop.txt
    if n < size(blockEvents,1)
        nextBlock = seq.getBlock(n+1);  
    else
        nextBlock = [];
    end

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
        %fprintf('n: %d, view: %d, sl: %d, echo: %d\n', n, view, sl, echo);

        adcCount = adcCount+1;
    end

    % create a TOPPE module struct from current Pulseq block
    moduleCandidate = pulsegeq.sub_block2module(block, n, systemGE, length(modules) + 1);

    % Is there an existing module that can be reused (scaled)?
    % Specifically, does one of the existing modules (elements of modules vector) have 
    % the same length waveform, the same non-empty rf/gx/gy/gz,
    % and the same value of 'hasRF' and 'hasADC', as moduleCandidate?
    isUnique = 1; 
    for p = 1:length(modules)
        if (modules(p).res == moduleCandidate.res ...
            & isempty(modules(p).rf) == isempty(moduleCandidate.rf) ...
            & isempty(modules(p).gx) == isempty(moduleCandidate.gx) ...
            & isempty(modules(p).gy) == isempty(moduleCandidate.gy) ...
            & isempty(modules(p).gz) == isempty(moduleCandidate.gz) ...
            & modules(p).hasRF  == moduleCandidate.hasRF ...
            & modules(p).hasADC == moduleCandidate.hasADC ...
            )
            isUnique = 0;
            break;   % break out of 'p' loop. 'p' now has the value of a module we'll reuse
        end
    end

    if isUnique
        % We found a unique block, so add it as a new module
        if arg.verbose
            fprintf('\tFound new module at block %d\n', n);
        end
        modules(end+1) = moduleCandidate;
        loop(n) = pulsegeq.sub_updateloopstruct([], block, nextBlock, systemGE, ...
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', length(modules));
        continue; % done, so move on to next block (n loop)
    end

    % moduleCandidate is not unique and has the same waveform length as an existing module, 
    % so now check to see if all waveform shapes in moduleCandidate match those in modules(ic).
    % If not, add waveform to module waveform group.

    tol = 1e-3;  % Shape is deemed equal if sum(abs(difference)) < tol
    ii = 1;
    isSameShape = [];
    for ax = {'rf','gx','gy','gz'};
        ax = cell2mat(ax);
        ch = block.(ax);
        if ~isempty(ch)
            for w = 1:modules(p).npulses  % Loop variable name 'w' here is consistent with usage in interpreter
                eval(sprintf('wav1 = modules(p).%s(:,w);', ax));
                eval(sprintf('wav2 = moduleCandidate.%s;', ax));
                isSameShape(ii,w) = norm(wav1-wav2,1) < tol;
            end
            ii = ii + 1;
        end
    end

    res = sum(isSameShape,1) == size(isSameShape,1);
    I = find(res==1);
    if ~isempty(I)
        % We found a set of RF/gradient waveforms in modularArr(p) with the same shapes as those in moduleCandidate,
        % so we'll reuse that and set 'mod' and 'wavnum' (waveform array column index) accordingly.
        wReuse = I(1);
        loop(n) = pulsegeq.sub_updateloopstruct([], block, nextBlock, systemGE, ...
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', p, 'wavnum', wReuse);
    else
        % Found a new set of shapes, so add this waveform set to modules(p)
        modules(p) = pulsegeq.sub_updatemodule(modules(p), block, n, systemGE);
        loop(n) = pulsegeq.sub_updateloopstruct([], block, nextBlock, systemGE, ... 
            'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', p, 'wavnum', modules(p).npulses);
    end

end

if arg.verbose
    fprintf(' done\n');
else
    fprintf('\n');
end


%% First, write each module to a .mod file.
%% Write modules.txt along the way.
if arg.verbose
    fprintf(1, 'Writing .mod files and modules.txt...\n');
end

% write modules.txt header
fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n', length(modules));
fprintf(fid,'wavfile_name    duration (us)     has_RF?     has_ADC?\n');

% loop through modules
for p = 1:length(modules)

    hasADC = modules(p).hasADC;
    hasRF  = modules(p).hasRF;

    if hasRF & hasADC
        error('Cannot transmit RF and acquire data in same block. Redesign the .seq file.');
    end

    % waveforms in modules are normalized, so now we scale to physical units
    rf = [];
    gx = [];
    gy = [];
    gz = [];
    if hasRF
        for w = 1:modules(p).npulses
            rfmax = 0;
            for n = 1:length(loop)
                if loop(n).mod == p & loop(n).wavnum == w
                    rfmax = max(loop(n).rfamp, rfmax);
                end
            end
            rf(:, w) = rfmax * modules(p).rf(:, w);
            RFmax(p, w) = rfmax;
            %modules(p).rf(:, w) = rfmax * modules(p).rf(:, w);
        end
    end

    for w = 1:modules(p).npulses
        gxmax = 0;
        gymax = 0;
        gzmax = 0;
        for n = 1:length(loop)
            if loop(n).mod == p & loop(n).wavnum == w
                gxmax = max(abs(loop(n).gxamp), gxmax);
                gymax = max(abs(loop(n).gyamp), gymax);
                gzmax = max(abs(loop(n).gzamp), gzmax);
            end
        end
        for ax = {'gx','gy','gz'}
            ax = cell2mat(ax);
            eval(sprintf('wav = modules(p).%s;', ax));
            if ~isempty(wav)
                eval(sprintf('%s(:,w) = %smax * modules(p).%s(:,w);', ax, ax, ax));
                %eval(sprintf('modules(p).%s(:,w) = %smax * modules(p).%s(:,w);', ax, ax, ax));
            end
        end
        GXmax(p, w) = gxmax;
        GYmax(p, w) = gymax;
        GZmax(p, w) = gzmax;
        %[p w rfmax gxmax gymax gzmax]
    end

    %fprintf('module %d: rfmax: %.3f, gxmax:%.2f, gymax:%.2f, gzmax:%.2f\n', p, rfmax, gxmax, gymax, gzmax);

    if arg.verbose
        fprintf('Creating .mod file number %d...\n', p);
    end

    % make sure waveforms start and end at zero, and are on a 4-sample boundary (toppe.writemod requires this)
    channels = {'rf','gx','gy','gz'};
    for ii=1:length(channels)
        eval(sprintf('wav = %s;', channels{ii}));
        if ~isempty(wav)
            [nt npulses] = size(wav);
            if any(wav(1,:) ~= 0)
                wav = [zeros(1,npulses); wav];
            end
            if any(wav(end,:) ~= 0)
                wav = [wav; zeros(1,npulses)];
            end
            [nt npulses] = size(wav);
            wav = toppe.makeGElength(wav);
        end
        eval(sprintf('%s = wav;', channels{ii}));
    end

    try
        warning('off');    % don't show message about padding waveforms
        toppe.writemod(systemGE, 'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', modules(p).ofname); 
        warning('on');
    catch ME
        error(sprintf('Error in writemod:\n%s', ME.message));
    end
        
    if arg.verbose
        fprintf('success\n');
    end

    % update entry in modules.txt

    fprintf(fid,'%s\t%d\t%d\t%d\t%d\n', modules(p).ofname, 0, hasRF, hasADC, modules(p).trigpos);    
end
fclose(fid);

if arg.verbose
    fprintf('done. Created %d .mod files.\n', p);
    toppe.plotmod('all');
end


%% Write scanloop.txt, which specifies the scan sequence (along with modules.txt and the .mod files).

% load .mod files
mods = toppe.tryread(@toppe.readmodulelistfile, 'modules.txt');

toppe.write2loop('setup', systemGE, 'version', arg.toppeVersion); 

for n = 1:length(loop)

    if isempty(loop(n).mod)
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

    % p and w follow naming convention in interpreter
    p = loop(n).mod;
    w = loop(n).wavnum;

    % RF scaling
    if modules(p).hasRF
        RFamplitude = loop(n).rfamp/RFmax(p,w);
    end

    % gradient scaling
    if GXmax(p,w) > 0
        Gamplitude(1) = loop(n).gxamp/GXmax(p,w);
    end
    if GYmax(p,w) > 0
        Gamplitude(2) = loop(n).gyamp/GYmax(p,w);
    end
    if GZmax(p,w) > 0
        Gamplitude(3) = loop(n).gzamp/GZmax(p,w);
    end

    RFphase  = loop(n).rfphs;
    DAQphase = loop(n).recphs;
    RFspoil  = false;
    RFoffset = loop(n).rffreq;    % Hz
    slice    = loop(n).slice;
    echo     = loop(n).echo + 1;  % write2loop starts indexing at 1
    view     = loop(n).view;
    view     = loop(n).view;
    Dabmodes = {'off','on'};
    dabmode  = Dabmodes{loop(n).dabmode+1};
    textra   = loop(n).textra*1e3;    % msec
    trigout  = loop(n).trigout;

    if textra < 0
        textraWarning = true;
        textra = 0;
    else
        textraWarning = false;
    end

    toppe.write2loop(modules(p).ofname, systemGE, ...
        'RFamplitude', RFamplitude, ...
        'Gamplitude',  Gamplitude, ...
        'slice',       slice, ...
        'echo',        echo, ...
        'view',        view, ...
        'dabmode',     dabmode, ...
        'RFphase',     RFphase, ...
        'DAQphase',    DAQphase, ...
        'textra',      textra, ...
        'RFoffset',    RFoffset, ...
        'waveform',    w, ...
        'trigout',     trigout);

end

if textraWarning
    fprintf(['\nWarning: requested textra < 0, which means that .seq sequence timing ', ...
        'is too tight to be directly converted to TOPPE --', ...
        ' ''textra'' set to zero in one or more scanloop.txt entries.\n']);
end

toppe.write2loop('finish', systemGE);

return

if arg.toppeVersion > 5
    % Write cores.txt, which defines the block groups
    blockGroups = [];
    for ie=1:length(loop)
        bgID = loop(ie).blockGroupID;
        modID = loop(ie).mod;
        if ~isempty(bgID)
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
for p = 1:length(modules)
    if modules(p).hasRF
        b1ScalingFile = modules(p).ofname;
    end
    if modules(p).hasADC
        readoutFile = modules(p).ofname;
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
for p = 1:length(modules)
    system(sprintf('tar rf %s %s', arg.tarFile, modules(p).ofname));
end

% clean up
system('rm toppeN.entry seqstamp.txt modules.txt scanloop.txt');
if arg.toppeVersion > 5
    system('rm cores.txt');
end
for p = 1:length(modules)
    system(sprintf('rm %s', modules(p).ofname));
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

