function seq = ge2seq(toppeTarFile, systemGE, systemSiemens, varargin)
% function seq = ge2seq(toppeTarFile, systemGE, systemSiemens, varargin)
%
% TOPPE to Pulseq file conversion.
%
% Inputs:
%  toppeTarFile       TOPPE .tar archive file containing the following TOPPE scan files:
%                     *.mod files       
%                     modules.txt     
%                     scanloop.txt    S
%                     If empty ([]), these files are asssumed to exist in the local path.
%  systemGE           struct specifying GE system specs, see toppe.systemspecs()
%  systemSiemens      struct specifying Siemens scanner hardware limits, see mr.opts()
% Options:
%  seqFile            Output .seq file name
%  FOV                [1 3] (m)
%  name               string
%  moduleListFile     Text file listing all .mod files. Default: 'modules.txt'.
%                     The .mod files listed must exist in the Matlab path.
%  loopFile           Text file specifying the MR scan loop. Default: 'scanloop.txt'
%  nt                 Only step through the first nt rows in scanloop.txt. Default: all rows.
%
% Examples:
%  >> pulsegeq.ge2seq('cal.tar', 'seqFile', 'cal.seq');
%  >> lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
%                    'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
%                    'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  
%  >> sys = systemspecs('maxSlew', 130, 'slewUnit', 'T/m/s');
%  >> ge2seq('cal.tar', sys, lims);
%

import pulsegeq.*


%% Parse inputs
% defaults
arg.seqFile        = 'out.seq';
arg.FOV            = [];
arg.name           = 'from_ge2seq';
arg.debug          = false;
arg.debugAdc       = false;
arg.moduleListFile = 'modules.txt';
arg.loopFile       = 'scanloop.txt';
arg.nt             = [];

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);


%% Untar files
if ~isempty(toppeTarFile)
    try
        system(sprintf('tar xf %s --wildcards *.mod modules.txt scanloop.txt', toppeTarFile));
    catch ME
        error(ME.message);
        return;
    end
else
    % The following TOPPE scan files must exist in current folder:
    % modules.txt, scanloop.txt, and the .mod files listed in modules.txt
end


%% Read TOPPE scan info
max_pg_iamp  = 2^15-2;    % max TOPPE/GE "instruction amplitude" (signed short int)
d      = toppe.tryread(@toppe.readloop,           arg.loopFile);         % scanloop array
modArr = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);   % module waveforms


%% set number of rows in scanloop.txt to step through
if isempty(arg.nt)
    nt = size(d,1);    % number of startseq calls
else
    nt = arg.nt;
end


%% Loop through scanloop.txt. Add each row as one Pulseq "block".
% Thankfully, no need to check for uniqueness since
% the Pulseq Matlab package already does that for us.

% initialize Pulseq sequence object
seq = mr.Sequence(systemSiemens);

raster = systemGE.raster;  % 4e-6 s

t = 0;  % running time 

for ii = 1:nt

    if ~mod(ii,100)
        for ib=1:60; fprintf('\b'); end;
        fprintf('ge2seq: parsing scan loop (%d of %d)', ii, nt);
    end

    module = modArr{d(ii,1)};

    % get number of discarded samples at beginning+end of RF waveform / ADC window
    nChop = [module.npre  module.res - module.rfres - module.npre];

    % get scaled waveforms (as row vectors)
    rfwav = (module.rf(:,1)).';  % full scale -- scaling done in makeArbitraryRf call
    gxwav = d(ii,4)/max_pg_iamp*(module.gx(:,1))';
    gywav = d(ii,5)/max_pg_iamp*(module.gy(:,1))';
    gzwav = d(ii,6)/max_pg_iamp*(module.gz(:,1))';

    % apply 3d rotation 
    Rv = d(ii,17:25)/max_pg_iamp;  % stored in row-major order
    R = reshape(Rv, 3, 3);
    G = R * [gxwav(:)'; gywav(:)'; gzwav(:)'];
    gxwav = G(1,:);
    gywav = G(2,:);
    gzwav = G(3,:);

    % convert to Pulseq units
    % rf:   Hz
    % grad: Hz/m
    %if ii==258; keyboard;  end;
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

    % additional delay at end of module
    textra = d(ii,14);   % us (microseconds)

    % Make Pulseq gradient structs (including all zero waveforms)
    if module.hasRF
        delay.grad = max(systemGE.start_core_rf*1e-6, systemSiemens.rfDeadTime);  % sec
    elseif module.hasDAQ
        delay.grad = max(systemGE.start_core_daq*1e-6, systemSiemens.adcDeadTime);  % sec
    else
        delay.grad = systemGE.start_core_grad*1e-6;
    end
    delay.grad = roundtoraster(delay.grad, systemSiemens.gradRasterTime);  % make multiple of 10us
    gx = mr.makeArbitraryGrad('x', gxwavPulseq, systemSiemens, 'delay', delay.grad);
    gy = mr.makeArbitraryGrad('y', gywavPulseq, systemSiemens, 'delay', delay.grad);
    gz = mr.makeArbitraryGrad('z', gzwavPulseq, systemSiemens, 'delay', delay.grad);

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
    strArg = getStrArg(hasg);    % 'gz' or 'gx,gy,gz' or... as appropriate

    if module.hasRF
        if nChop(2) < ceil(systemSiemens.rfRingdownTime/raster)
            error(sprintf('%s: RF ringdown occurs past end of gradient -- increase nChop(2)', module.fname));
            %nChop(2) = 48;
            %warning(sprintf('%s: RF ringdown occurs past end of gradient -- nChop(2) set to 48 samples', module.fname));
        end

        phaseOffset = d(ii,12)/max_pg_iamp*pi;  % radians
        freqOffset  = d(ii,15);                 % Hz

        % get nominal/full flip angle from .mod file header, and flip angle scaling for this block
        nominalFlip = module.paramsfloat(16)/180*pi;   %  assumes that flip angle is stored in .mod file header
        flipScaling = d(ii, 2)/max_pg_iamp;

        % start of RF waveform. Don't include psd_rf_wait.
        delay.rf = delay.grad + nChop(1)*raster;
        %delay.rf = max(0, systemGE.start_core_rf*1e-6 - nChop(1)*raster) + ...
        %    systemGE.myrfdel*1e-6 + nChop(1)*raster ;
        %delay.rf = roundtoraster(delay.rf, systemSiemens.gradRasterTime); 

        % rf object
        iStart = 1 + round(nChop(1)*raster / systemSiemens.rfRasterTime);
        iStop = iStart - 1 + round( (length(rfwav) - sum(nChop))*raster / systemSiemens.rfRasterTime);
        rf = mr.makeArbitraryRf(rfwavPulseq(iStart:iStop), flipScaling*nominalFlip, ...
            'PhaseOffset', phaseOffset, ...
            'FreqOffset', freqOffset, ...
            'system', systemSiemens, ...
            'delay', delay.rf);

        % create block. Include gradient objects if any.
        if isempty(strArg)
            seq.addBlock(rf); %, adcPad);
        else
            %eval( sprintf( 'seq.addBlock(rf, %s, adcPad)', strArg) ); % TODO: remove adcPad
            eval( sprintf( 'seq.addBlock(rf, %s)', strArg) ); 
        end

        % Delay after end of RF waveform
        % Don't add psd_rf_wait (=myrfdel) here.
        %postDelay = max(0, systemGE.myrfdel*1e-6 - raster*nChop(2)) ...  % "coredel" in toppe.plotseq()
        %    + (systemGE.timetrwait + systemGE.tminwait + systemGE.timessi + textra)*1e-6;
        postDelay = (systemGE.timetrwait + systemGE.tminwait + systemGE.timessi + textra)*1e-6;

        if arg.debug
            clf;
            subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal)))); ylabel('Hz');
            subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal)))); ylabel('rad');
        end
    elseif module.hasDAQ
        % Acquire with same raster/dwell time as TOPPE (4us)
        % Make number of samples a multiple of 10 to enforce
        % an even number of samples and duration on 10us boundary.
        ntmp = numel(gxwav) - sum(nChop); 
        nADC = ntmp - mod(ntmp, 10);  % number of data samples to acquire

        phaseOffset = d(ii,13)/max_pg_iamp*pi;          % radians

        % start of ADC window 
        delay.adc = delay.grad + nChop(1)*raster;
        %delay.adc = max(0, delay.grad - nChop(1)*raster) ...
        %    + systemGE.daqdel*1e-6 + nChop(1)*raster;
        %delay.adc = roundtoraster(delay.adc, systemSiemens.gradRasterTime); 

        if delay.adc < systemSiemens.adcDeadTime
            warning(sprintf('Requested ADC start time < systemSiemens.adcDeadTime. Extended to %.3e.', systemSiemens.adcDeadTime));
            delay.adc = systemSiemens.adcDeadTime;
        end

        % ADC object
        adc = mr.makeAdc(nADC, ...
            'Dwell', raster, ...
            'PhaseOffset', phaseOffset, ...
            'system', systemSiemens, ...
            'delay', delay.adc);

        % create block. Include gradient objects if any.
        if isempty(strArg)
            seq.addBlock(adc);
        else
            eval( sprintf( 'seq.addBlock(%s, adc)', strArg) );
        end

        % Delay after end of ADC window 
        % Don't add daqdel.
        postDelay = (systemGE.timetrwait + systemGE.tminwait + systemGE.timessi + textra)*1e-6;
        %postDelay = max(0, systemGE.daqdel*1e-6 - raster*nChop(2)) ...
        %    + (systemGE.timetrwait + systemGE.tminwait + systemGE.timessi + textra)*1e-6;

        % Round block duration to 10us boundary (this seems to be a Pulseq requirement)
        %blk = seq.getBlock(blkIndex);
        %.blockDuration = roundtoraster(adc.blockDuration);
    else
        % Create block containing only gradients (if any)
        if ~isempty(strArg)
            eval( sprintf( 'seq.addBlock(%s)', strArg) );
        end
        postDelay = (systemGE.timetrwait + systemGE.tminwait + systemGE.timessi + textra)*1e-6;
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

    % Add delay at end of block.
    % Delay must (?) be in multiples of gradient raster times.
    postDelay = roundtoraster(postDelay, systemSiemens.gradRasterTime); 
    del = mr.makeDelay(postDelay);
    seq.addBlock(del);
    clear postDelay; 
end
fprintf('\n');


%% Check sequence timing and write to file
fprintf('Checking Pulseq timing... ');
[ok, error_report]=seq.checkTiming;
if (ok)
    if ~isempty(arg.FOV)
        seq.setDefinition('FOV', arg.FOV);
    end
    seq.setDefinition('Name', arg.name);
    seq.write(arg.seqFile);
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%rep = seq.testReport;
%fprintf([rep{:}]); 

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
