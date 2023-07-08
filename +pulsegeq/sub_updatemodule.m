function module = sub_updatemodule(module, block, blockid, system)
% Add waveforms from 'block' to an existing module waveform array
%
% Add the (normalized) waveforms contained in the block number 'blockid' 
% as the last column of module.rf/module.<grad>.
%
% Inputs:
%   module           see sub_block2module()
%   block            Pulseq block
%   blockid          Pulseq block is seq.getBlock(blockid)
%   system           TOPPE system struct

import pulsegeq.*

raster  = system.raster;   % sec

if isempty(module.blockids)
    module.blockids(1) = blockid;
else
    module.blockids(end+1) = blockid;
end

if ~isempty(block.adc) & ~isempty(block.rf)
    error('Block/module can not be both RF transmit and receive');
end

module.npre = 0;   % RF/ADC delay (number of 4us samples)
module.rfres = 0;  % Duration of RF/ADC window (number of 4us samples)

% RF
if ~isempty(block.rf)
    module.hasRF = 1;

    % interpolate to GE raster time (4us)
    %rf = downsample(block.rf.signal, round(raster/1e-6));     % downsample from 1us to 4us (GE raster time)
    tge = block.rf.t(1) : system.raster : block.rf.t(end);
    rf = interp1(block.rf.t, block.rf.signal, tge, 'linear', 'extrap');     % downsample from 1us to 4us (GE raster time)

    % make length even
    %rf = toppe.makeGElength(rf(:));

    module.rfres = length(rf);

    % pad with zeros during delay (the .mod file format requires it; interpreter ignores these zeros)
    module.npre = round(block.rf.delay/raster);
    rf = [linspace(0, 0, module.npre)'; rf(:)];

    % Normalize and add to waveforms
    rf = rf/max(abs(rf(:)));
    module.rf = sub_addwav(module.rf, rf);
end

% Gradients
for ax = {'gx','gy','gz'};
    ax = ax{1};
    grad = block.(ax);
    if ~isempty(grad)
        if strcmp(grad.type, 'grad')    % arbitrary shape
            % interpolate to GE raster time (4us)
            tge = raster/2:raster:grad.tt(end);
            wav = interp1(grad.tt, grad.waveform, tge, 'linear', 'extrap'); 
            wav = wav/100/system.gamma;     % Gauss/cm
        else
            % trapezoid
            wav = sub_trap2shape(grad, system.raster);     % Gauss/cm
        end

        % check slew
        peakSlew = max(max(diff(wav,1)/(raster*1e3),[],1));
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
    module.npre = round(block.adc.delay/raster);

    % ADC duration = rfres*4us
    adcDuration = block.adc.dwell*block.adc.numSamples;
    module.rfres = 2*ceil(adcDuration/system.raster/2);
end

% If ADC block without gradients, create 'dummy' waveform to keep writemod happy
% (needs at least one non-zero waveform)
if module.hasADC & length([module.rf(:); module.gx(:); module.gy(:); module.gz(:)]) == 0
    module.gx = 0.01*ones(round((block.adc.delay + block.adc.dwell*block.adc.numSamples)/raster), 1);
end

% total number of 4us samples in module
module.res = max([ length(module.rf) length(module.gx) ...
           length(module.gy) length(module.gz) ...
           module.rfres]);

if module.rfres == 0
    module.rfres = module.res;
end

% pad with zeros as needed to ensure equal length of all (non-empty) waveforms
npulses = 0;
if ~isempty(module.rf)
    wav = module.rf;
    npulses = size(wav,2);
    module.rf = [wav; zeros(module.res-size(wav,1), size(wav,2))];
end
gradChannels = {'gx','gy','gz'};
for ii=1:length(gradChannels)
    eval(sprintf('wav = module.%s;', gradChannels{ii}));
    if ~isempty(wav)
        npulses = size(wav,2);
        wav = [wav; zeros(module.res-size(wav,1), size(wav,2))];
        eval(sprintf('module.%s= wav;', gradChannels{ii}));
    end
end

module.npulses = npulses;

return
