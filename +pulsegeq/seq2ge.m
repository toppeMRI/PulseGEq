function [moduleArr loopStructArr] = seq2ge(seqarg, varargin)
% function seq2ge(seqarg, varargin)
%
% Convert a Pulseq file (http://pulseq.github.io/) to a set of TOPPE files
% that can be executed on GE MR scanners. 
%
% See https://toppemri.github.io/ for more info on TOPPE.
%
% This script writes the following files to disk:
%   *.mod:            One .mod file corresponds to one "unique" block (see below)
%   modules.txt       List of .mod files, and flags indicating whether each .mod file corresponds to an RF/ADC/(gradients only) module
%   scanloop.txt      Sequence of instructions for the entire scan (waveform amplitudes, ADC instructions, etc)
%
% Inputs:
%   seqarg            Either a Pulseq file name, or an mr.Sequence object.
% Input options:
%   system            struct        Contains GE and TOPPE system specs, including TOPPE version. See +toppe/systemspecs.m
%   toppeVersion      string        'v2' (default) or 'v3'
%   verbose           boolean       Default: false
%   debug             boolean       Display detailed info about progress (default: false)
%   pulseqVersion     string        'v1.3.0' (default) or 'v1.2.1'
%   tarFile           string        default: 'toppeScanFiles.tar'
%   blockStop         int           end at this block in the .seq file (for testing)
%
% Usage examples:
%   >> seq2ge('../examples/2DFLASH.seq', 'toppeVersion', 'v3');
%   >> seq2ge('../examples/2DFLASH_v1.2.1.seq', 'pulseqVersion', 'v1.2.1');
%
%   >> system = toppe.systemspecs('maxSlew',200,'slewUnit','T/m/s','maxGrad',50','gradUnit','mT/m');
%   >> seq2ge('2DFLASH.seq', 'system', system, 'verbose', true);
%
%   >> seq = mr.Sequence();
%   >> seq.read('2DFLASH.seq');
%   >> seq2ge(seq, 'system', system);
%

%import pulsegeq.*

%% parse inputs
% Defaults
arg.system  = toppe.systemspecs();
arg.toppeVersion = 'v3';
arg.verbose = false;
arg.debug = false;
arg.pulseqVersion = 'v1.3.1';
arg.tarFile = 'toppeScanFiles.tar';
arg.blockStop = [];
arg.ibstart = 1;

%  systemSiemens      struct containing Siemens system specs. 
%                        .rfRingdownTime     Default: 30e-6   (sec)
%                        .rfDeadTime         Default: 100e-6  (sec)
%                        .adcDeadTime        Default: 20e-6   (sec)
%arg.systemSiemens = struct('rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6,'adcDeadTime', 20e-6);

% Substitute specified system values as appropriate (from MIRT toolbox)
arg = toppe.utils.vararg_pair(arg, varargin);

switch arg.toppeVersion
	case 'v2' 
		nCols = 16;   % number of columns in scanloop.txt
	case 'v3' 
		nCols = 25;   % number of columns in scanloop.txt
	otherwise
		error('Please use TOPPE v2 or v3');
end

switch arg.pulseqVersion
	case 'v1.2.1'
		nEvents = 6;   % number of events per block (number of columns in .seq file)
	case 'v1.3.0'
		nEvents = 7;
	case 'v1.3.1'
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

%% Loop through blocks and build 'moduleArr' and 'loopStructArr'

% 'moduleArr' struct array
% Find blocks that are unique in terms of waveforms and timing (i.e., waveform amplitudes, RF/ADC phase, etc can differ),
% and fill 'moduleArr' struct array accordingly. 
% Each entry of 'moduleArr' is a struct containing all waveforms belonging to one module (.mod file), and other module info.
% The usage of the word "module" here is consistent with its usage in TOPPE.
% For now, ignore 'EXT' blocks (TODO)

% 'loopStructArr' struct array
% Each entry in this array contains information needed to fill out one row of scanloop.txt.

if arg.verbose
	fprintf('Filling moduleArr struct array, and loopStructArr array.\n' );
end

% get contents of [BLOCKS] section
blockEvents = cell2mat(seq.blockEvents);
blockEvents = reshape(blockEvents, [nEvents, length(seq.blockEvents)]).'; % hardcoded for as long as Pulseq does not include another element

if ~isempty(arg.blockStop)
	blockEvents = blockEvents(1:arg.blockStop, :);
end

% First entry in 'moduleArr' struct array
block = seq.getBlock(arg.ibstart);
if ~isempty(block.delay)
	error('First block can''t contain a delay. Edit the .seq file.');
end
moduleArr(1) = pulsegeq.sub_block2module(block, arg.ibstart, arg.system, 1);

% First entry in 'loopStructArr' struct array (first block is by definition a module)
nextblock = seq.getBlock(arg.ibstart+1);   % needed to set 'textra' in scanloop.txt
loopStructArr(1) = pulsegeq.sub_updateloopstruct([], block, nextblock, arg.system, 'mod', 1);

% data frames (in Pfile) are stored using indeces 'slice', 'echo', and 'view' 
sl = 1;
view = 1;
echo = 0; 
adcCount = 0;

% h = waitbar(0,'Looping through blocks and looking for uniqueness...');
for ib = (arg.ibstart+1):size(blockEvents,1)
	if ~mod(ib, 500)
	%	waitbar(ib/size(blockEvents,1),h)
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
	if isfield(nextblock, 'trig') 
		nextblock = [];
	end

	%if ~isempty(block.delay) & isempty(block.rf) & isempty(block.adc) ...
	if isempty(block.rf) & isempty(block.adc) ...
		& isempty(block.gx) & isempty(block.gy) & isempty(block.gz) ...
		| isfield(block, 'trig')  % ignore trigger (ext) blocks for now. TODO
		continue;  % Done, move on to next block 
		% Pure delay blocks are accounted for in 'textra' in the previous row in scanloop.txt)
	end

	% set slice/echo/view indeces (if block is an acquisition block)
	% view = 1, ..., system.maxView
	% sl   = 1, ..., system.maxSlice
	if ~isempty(block.adc)
		view = mod(adcCount, arg.system.maxView) + 1;
		sl   = floor(adcCount/arg.system.maxView) + 1;
		if sl > arg.system.maxSlice;
			error(sprintf('max number of slices ecxeeded (%d)', arg.system.maxSlice));
		end
		echo = floor(adcCount/(arg.system.maxView*arg.system.maxSlice));
		if echo > arg.system.maxEcho
			error(sprintf('max number of echoes ecxeeded (%d)', arg.system.maxEcho));
		end
		%fprintf('ib: %d, view: %d, sl: %d, echo: %d\n', ib, view, sl, echo);

		adcCount = adcCount+1;
	end

	% create a TOPPE module struct from current Pulseq block
	modCandidate = pulsegeq.sub_block2module(block, ib, arg.system, length(moduleArr) + 1);

	% Is there an existing module that can be 'reused'?
	% Specifically, does one of the existing modules (elements of moduleArr) have the same length waveform, 
	% the same non-empty rf/gx/gy/gz, and the same value of 'hasADC', as modCandidate?
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
		loopStructArr(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, arg.system, ...
			'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', length(moduleArr));
		continue; % done, so move on to next block
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
		% so we'll 'reuse' that and set 'mod' and 'wavnum' (waveform array column index) accordingly.
		iWavReuse = I(1);
		loopStructArr(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, arg.system, ...
			'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic, 'wavnum', iWavReuse);
	else
		% Found a new set of shapes, so add this waveform set to moduleArr(ic)
		moduleArr(ic) = pulsegeq.sub_updatemodule(moduleArr(ic), block, ib, arg.system);
		loopStructArr(ib) = pulsegeq.sub_updateloopstruct([], block, nextblock, arg.system, ... 
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
		toppe.writemod('system', arg.system, 'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', moduleArr(ic).ofname); 
		warning('on');
	catch ME
		error(sprintf('Error in writemod:\n%s', ME.message));
	end
		
	if arg.verbose
		fprintf('success\n');
	end

	% update entry in modules.txt
	fprintf(fid,'%s\t%d\t%d\t%d\n', moduleArr(ic).ofname, 0, hasrf, hasadc);	
end
fclose(fid);

if arg.verbose
	fprintf('done. Created %d .mod files.\n', ic);
	toppe.plotmod('all');
end


%% Write scanloop.txt, which specifices the scan sequence (along with modules.txt and the .mod files).

% load .mod files
mods = toppe.utils.tryread(@toppe.readmodulelistfile, 'modules.txt');

toppe.write2loop('setup', 'version', str2num(arg.toppeVersion(2))); 

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
	toppe.write2loop(moduleArr(iMod).ofname, ...
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

toppe.write2loop('finish');

if arg.verbose
	fprintf(' done\n');
end

%% Put TOPPE files in a .tar file (for convenience)
system(sprintf('tar cf %s modules.txt scanloop.txt', arg.tarFile));
for ic = 1:length(moduleArr)
	system(sprintf('tar rf %s %s', arg.tarFile, moduleArr(ic).ofname));
end

return;

% list archive contents
if arg.verbose
	fprintf('\nCreated %s containing the following files:\n', arg.tarFile);
	system(sprintf('tar tf %s', arg.tarFile));
end

% clean up
system('rm modules.txt scanloop.txt');
for ic = 1:length(moduleArr)
	system(sprintf('rm %s', moduleArr(ic).ofname));
end

if ~arg.verbose
	fprintf('\n');
end

fprintf('Remember to rename one of the .mod files to ''tipdown.mod'', and another to ''readout.mod''\n');

return

%% End of main script

