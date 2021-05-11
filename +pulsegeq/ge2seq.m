function seq = ge2seq(toppeTarFile, varargin)
% function seq = ge2seq(toppeTarFile, varargin)
%
% TOPPE to Pulseq file conversion.
% Inputs:
%  toppeTarFile       TOPPE .tgz archive file containing the following:
%                     *.mod:            One .modfile corresponds to one "unique" block (see below)
%                     modules.txt       List of .mod files, and flags indicating whether the .wav file is RF/ADC/(gradients only)
%                     scanloop.txt      Sequence of instructions for the entire scan (waveform amplitudes, ADC instructions, etc)
% Options:
%  seqFile            Output .seq file name
%  moduleListFile     Text file listing all .mod files. Default: 'modules.txt'.
%                     The .mod files listed must exist in the Matlab path.
%  loopFile           Text file specifying the MR scan loop. Default: 'scanloop.txt'
%  system             struct specifying Siemens scanner hardware limits, see mr.opts
%  systemGE           struct specifying GE system specs, see +toppe/systemspecs.m
%
% Examples:
%  >> pulsegeq.ge2seq('cal.tar', 'seqFile', 'cal.seq');
%  >> lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
%                    'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
%                    'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  
%  >> sys = systemspecs('maxSlew', 130, 'slewUnit', 'T/m/s');
%  >> ge2seq('cal.tar', 'system', lims, 'systemGE', sys);
%

import pulsegeq.*

%% Parse inputs and set system values
% defaults
arg.seqFile        = 'out.seq';
arg.debug          = false;
arg.debugAdc       = false;
arg.moduleListFile = 'modules.txt';
arg.loopFile       = 'scanloop.txt';
arg.systemGE = toppe.systemspecs('addDelays', false);   % don't add delays before creating Pulseq blocks

raster = arg.systemGE.raster;  % 4e-6 s. For RF, gradients, and ADC.

arg.system = mr.opts('rfRasterTime', 1e-6, 'gradRasterTime', 10e-6, ...
                     'rfDeadTime', 100e-6, 'rfRingdownTime', 20e-6, ...
                     'adcDeadTime', 20e-6);  

%arg.system = mr.opts('maxGrad', 40, 'GradUnit', 'mT/m',...
%                     'maxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
%                     'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

% system struct to be used when creating Pulseq blocks
lims = arg.system;

% Define delays to pass to 'plotseq' call
%rfdel = max(arg.system.toppe.myrfdel   = round(arg.systemSiemens.rfDeadTime*1e6);       % RF delay (us)
if 0
arg.system.toppe.daqdel     = max(round(arg.systemSiemens.adcDeadTime*1e6 - arg.system.toppe.daqdel),    0);      % ADC delay (us)
arg.system.toppe.timetrwait = max(round(arg.systemSiemens.rfRingdownTime*1e6 - arg.system.toppe.daqdel), 0);   % (us)
arg.system.toppe.myrfdel    = round(arg.systemSiemens.rfDeadTime*1e6);       % RF delay (us)
arg.system.toppe.daqdel     = round(arg.systemSiemens.adcDeadTime*1e6);      % ADC delay (us)
arg.system.toppe.timetrwait = round(arg.systemSiemens.rfRingdownTime*1e6);   % (us)
end

% Untar files
try
	system(sprintf('tar xf %s', toppeTarFile));
catch ME
	error(ME.message);
	return;
end

% Read TOPPE scan info
max_pg_iamp  = 2^15-2;                  % max TOPPE/GE "instruction amplitude" (signed short int)
d      = toppe.utils.tryread(@toppe.readloop,           arg.loopFile);         % scanloop array
modArr = toppe.utils.tryread(@toppe.readmodulelistfile, arg.moduleListFile);   % module waveforms

% clean up TOPPE files
if false
system('rm modules.txt scanloop.txt');
for ic = 1:length(modArr)
	system(sprintf('rm %s', modArr{ic}.fname));
end
end


%% Loop through scanloop.txt. Add each row as one Pulseq "block".

% initialize Pulseq sequence object
seq = mr.Sequence(lims);

nt = size(d,1);    % number of startseq calls
for ii = 1:nt

	if ~mod(ii,100)
		for ib=1:60; fprintf('\b'); end;
		fprintf('Parsing scan loop: %d of %d', ii, nt);
	end

	module = modArr{d(ii,1)};

	% get waveforms and delay for one row (one startseq call)
	% rf: Gauss; gradients: Gauss/cm; tdelay: microsec
	[~, ~, ~, ~, rfwav, gxwav, gywav, gzwav, tdelay] = toppe.plotseq(ii, ii, ...
		'loopArr', d, 'mods', modArr, 'doDisplay', false, 'system', arg.systemGE);  

	% pulseq likes row vectors
	rfwav = rfwav(:)';
	gxwav = gxwav(:)';
	gywav = gywav(:)';
	gzwav = gzwav(:)';

	% convert to Pulseq units
	% rf:   Hz
   % grad: Hz/m
	rfwavPulseq = rf2pulseq(rfwav,raster,seq);
	gxwavPulseq = g2pulseq( gxwav,raster,seq);
	gywavPulseq = g2pulseq( gywav,raster,seq);
	gzwavPulseq = g2pulseq( gzwav,raster,seq);

	% ensure equal duration (interpolation to Pulseq rastertimes can result in unequal duration)
	trf   = length(rfwavPulseq) * seq.rfRasterTime;
	tgrad = length(gxwavPulseq) * seq.gradRasterTime;
	ngradextra = ceil((trf-tgrad)/seq.gradRasterTime);
	gxwavPulseq = toppe.utils.makeevenlength( [gxwavPulseq zeros(1, ngradextra)] );
	gywavPulseq = toppe.utils.makeevenlength( [gywavPulseq zeros(1, ngradextra)] );
	gzwavPulseq = toppe.utils.makeevenlength( [gzwavPulseq zeros(1, ngradextra)] );
	tgrad  = length(gxwavPulseq) * seq.gradRasterTime;
	rfwavPulseq = [rfwavPulseq zeros(1,round((tgrad-trf)/seq.rfRasterTime))];

	% Make Pulseq gradient structs (even all zero waveforms)
	gx = mr.makeArbitraryGrad('x', gxwavPulseq, lims);
	gy = mr.makeArbitraryGrad('y', gywavPulseq, lims);
	gz = mr.makeArbitraryGrad('z', gzwavPulseq, lims);

	% bitmask indicating non-zero gradients
	hasg = 0;   
	if ~all(gxwavPulseq == 0)
		hasg = bitset(hasg,1);
	end
	if ~all(gywavPulseq == 0)
		hasg = bitset(hasg,2);
	end
	if ~all(gzwavPulseq == 0)
		hasg = bitset(hasg,3);
	end
	strArg = getStrArg(hasg);        % 'gz' or 'gx,gy,gz' or... as appropriate

	freqOffset  = d(ii,15);                         % Hz

	if module.hasRF
		phaseOffset = d(ii,12)/max_pg_iamp*pi;          % radians
		flip = pi; % d(ii,2)/max_pg_iamp*pi;

		rf = mr.makeArbitraryRf(rfwavPulseq, flip, 'FreqOffset', freqOffset, ...
			'PhaseOffset', phaseOffset, 'system', lims);

		gx.delay = lims.rfDeadTime; % - 10e-6;   %
		gy.delay = lims.rfDeadTime;
		gz.delay = lims.rfDeadTime;

		if isempty(strArg)
			seq.addBlock(rf);
		else
			eval( sprintf( 'seq.addBlock(rf, %s)', strArg) );
		end

		if arg.debug
			clf;
			subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal)))); ylabel('Hz');
			subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal)))); ylabel('rad');
		end
	elseif module.hasDAQ
		phaseOffset = d(ii,13)/max_pg_iamp*pi;          % radians
		if arg.debugAdc
			% delay and shorten adc window
			nadc = 2*round(numel(gx.waveform)/4);
			delay = round(nadc/4)*seq.gradRasterTime;
			adc = mr.makeAdc(nadc, lims, 'Dwell', raster, 'delay', delay, ...
				'freqOffset', freqOffset, 'phaseOffset', phaseOffset);
		else
			adc = mr.makeAdc(numel(gxwav), lims, 'Dwell', raster, 'delay', 0,...
				'freqOffset', freqOffset, 'phaseOffset', phaseOffset);
		end

		gx.delay = lims.adcDeadTime;
		gy.delay = lims.adcDeadTime;
		gz.delay = lims.adcDeadTime;

		if isempty(strArg)
			seq.addBlock(adc);
		else
			eval( sprintf( 'seq.addBlock(%s, adc)', strArg) );
		end
	else
		if ~isempty(strArg)
			eval( sprintf( 'seq.addBlock(%s)', strArg) );
		end
	end

	if arg.debug
		if ~module.hasRF
			clf;
		end
		subplot(2,2,[3 4]); 
		hold on
		if bitget(hasg, 1)
			plot(gx.waveform,'r'); ylabel('Hz/m'); hold on; 
		end
		if bitget(hasg, 2)
			plot(gy.waveform,'g'); hold on;
		end
		if bitget(hasg, 3)
			plot(gz.waveform,'b'); 
		end
		if ~hasg
			clf(sfh);
		end
			
		hold off; %title(sprintf('max = %f', max([gx.waveform gy.waveform gz.waveform])));
		input('press any key to continue');
	end

	% add delay block
	if tdelay > 12;    % minimum duration of wait pulse in TOPPE
		del = mr.makeDelay(roundtoraster(tdelay*1e-6, raster)); % delay also needs to be in multiples of raster times
		seq.addBlock(del);
	end
end
fprintf('\n');


%% Check whether the timing of the sequence is correct
fprintf('Checking Pulseq timing... ');
[ok, error_report]=seq.checkTiming;

if (ok)
	fprintf('Timing check passed successfully\n');
	seq.write(arg.seqFile);
else
	fprintf('Timing check failed! Error listing follows:\n');
	fprintf([error_report{:}]);
	fprintf('\n');
end

return;



%% helper function: get gradient arguments (as string) to pass to seq.addBlock()
function argStr = getStrArg(hasg)

switch hasg
	case 0
		argStr = '';  % no gradients
	case 1
		argStr = 'gx'; 
	case 2
		argStr = 'gy'; 
	case 4
		argStr = 'gz'; 
	case 3
		argStr = 'gx, gy'; 
	case 5
		argStr = 'gx, gz'; 
	case 6
		argStr = 'gy, gz'; 
	case 7
		argStr = 'gx, gy, gz'; 
end

return;
