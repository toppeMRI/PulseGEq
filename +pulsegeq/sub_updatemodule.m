function module = sub_updatemodule(module, block, blockid, system)
% Add waveforms from 'block' to an existing module waveform array
%
% Add the (normalized) waveforms contained in the block number 'blockid' 
% as the last column of module.rf/module.<grad>.
%
% Inputs:
%   module           see sub_block2module()
%   blockid          Pulseq block is seq.getBlock(blockid)
%   system           TOPPE system struct

import pulsegeq.*

dt  = system.raster;   % sec

if isempty(module.blockids)
    module.blockids(1) = blockid;
else
    module.blockids(end+1) = blockid;
end

if ~isempty(block.adc) & ~isempty(block.rf)
    error('Block can not be both RF transmit and receive');
end

% RF
if ~isempty(block.rf)
    module.hasRF = 1;

    % interpolate to GE raster time (4us)
    %rf = downsample(block.rf.signal, round(dt/1e-6));     % downsample from 1us to 4us (GE raster time)
    tge = block.rf.t(1) : system.raster : block.rf.t(end);
    rf = interp1(block.rf.t, block.rf.signal, tge, 'linear', 'extrap');     % downsample from 1us to 4us (GE raster time)
    %if ( length(block.rf.signal) > 24 )
    %else
    %   rf = block.rf.signal;                               % assumed to be already decimated in terms of dt_ge
    %   warning('rf waveform is < 24 points and will not be interpolated to GE raster time');
    %end

    % add delay (pad with zeros)
    rf = [linspace(0, 0, round(block.rf.delay/dt))'; rf.'];

    % Normalize shapes and add to waveform group
    rf = rf/max(abs(rf(:)));
    module.rf = sub_addwav(module.rf, rf);    % Normalized to amplitude 1
end

% Gradients
for ax = {'gx','gy','gz'};
    ax = ax{1};
    grad = block.(ax);
    if ~isempty(grad)
        if strcmp(grad.type, 'grad')    % arbitrary shape
            % must start and end on zero
            grad.raster = grad.tt(2) - grad.tt(1);
            if abs(grad.waveform(1)) > 0
                grad.waveform = [0; grad.waveform];
                grad.tt = [grad.raster/2; (grad.tt + grad.raster)];
            end
            if abs(grad.waveform(end)) > 0
                % Add 2 sample ramp to zero.
                % This is a hack to deal with pathological case in some .seq files.
                grad.waveform = [grad.waveform; grad.waveform(end)/2; 0];
                grad.tt = [grad.tt; [grad.tt(end) + grad.raster*[1;2]]];
            end

            % interpolate to GE raster time
            % TODO: shift by dt/2?
            tge = 0:dt:(max(grad.tt)-0*dt);
            wav = interp1(grad.tt, grad.waveform, tge, 'linear', 'extrap');   % interpolate to GE raster time (4us)
            wav(isnan(wav)) = 0;                         % must be due to interp1
            if wav(1) > 0
                wav = [0; wav(:)];
            end
            if wav(end) > 0
                wav = [wav(:); 0];
            end
            wav = wav/100/system.gamma;                  % Gauss/cm
        else
            % trapezoid
            wav = sub_trap2shape(grad, system.raster);     % Gauss/cm
        end

        % check slew
        peakSlew = max(max(diff(wav,1)/(dt*1e3),[],1));
        if peakSlew > system.maxSlew
            %[peakSlew system.maxSlew blockid module.modnum size(module.(ax),2)+1]
            error(sprintf('slew rate violation at blockid %d (%s) (%.1f%%)', blockid, ax, peakSlew/system.maxSlew*100));
        end 

        % normalize and add to module 
        if norm(wav) > 0
            wav = wav/max(abs(wav(:)));
        end
        module.(ax) = sub_addwav(module.(ax), wav(:));
    end
end

% ADC
if ~isempty(block.adc)
    module.hasADC = 1;
    %module.ofname = 'readout.mod';
    tend = block.adc.delay + block.adc.dwell*block.adc.numSamples;  % sec
    nAdc = round(tend/dt);
else
    nAdc = 0;
end

% if ADC block without gradients, create 'dummy' waveform to keep writemod happy
if module.hasADC & length([module.rf(:); module.gx(:); module.gy(:); module.gz(:)]) == 0
    module.gx = 0.01*ones(round((block.adc.delay + block.adc.dwell*block.adc.numSamples)/dt), 1);
end

% store waveform length (useful for comparing blocks)
res = max([ length(module.rf) length(module.gx) ...
           length(module.gy) length(module.gz) ...
           nAdc]);

% pad with zeros as needed to ensure equal length of all (non-empty) waveforms
npulses = 0;
if ~isempty(module.rf)
    wav = module.rf;
    npulses = size(wav,2);
    module.rf = [wav; zeros(res-size(wav,1), size(wav,2))];
end
gradChannels = {'gx','gy','gz'};
for ii=1:length(gradChannels)
    eval(sprintf('wav = module.%s;', gradChannels{ii}));
    if ~isempty(wav)
        npulses = size(wav,2);
        wav = [wav; zeros(res-size(wav,1), size(wav,2))];
        eval(sprintf('module.%s= wav;', gradChannels{ii}));
    end
end

module.res = res;
module.npulses = npulses;

% Set npre (RF/ADC delay, in number of 4us samples),
% and rfres (duration of RF/ADC, in number of 4us samples)
if module.hasADC
    module.npre = round(block.adc.delay/dt);
    module.rfres = nAdc;
elseif module.hasRF
    module.npre = round(block.rf.delay/dt);
    module.rfres = length(module.rf);
else
    module.npre = 0;
    module.rfres = module.res;
end

% trigger out
if isfield(block, 'trig')
    module.trigpos = round(block.trig.delay*1e6);   % us
else
    module.trigpos = -1;    % no trigger
end

return
