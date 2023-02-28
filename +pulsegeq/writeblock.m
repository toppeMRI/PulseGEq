function writeblock(fid, blkID, blk, sys)
%
% Write a Pulseq block to a custom machine-readable file format 
% for the PulseGEq interpreter.

sys.gamma = 4.2576e3;  % Hz/Gauss

% order is important
fwrite(fid, blkID, 'int16');
fwrite(fid, round(blk.blockDuration*1e6), 'int16');
sub_writerf(fid, blk.rf, sys);
sub_writegrad(fid, blk.gx, sys);
sub_writegrad(fid, blk.gy, sys);
sub_writegrad(fid, blk.gz, sys);
sub_writeadc(fid, blk.adc, sys);

return


function sub_writerf(fid, rf, sys)

C = pulsegeq.constants;

% type
if isempty(rf)
    fwrite(fid, C.NULL, 'int16');
    return
end
fwrite(fid, C.ARBITRARY, 'int16');

% downsample from 1us to 4us raster
rf.signal = rf.signal(2:4:end);

% delay
fwrite(fid, round(rf.delay*1e6), 'int16');     % us

% amplitude
amp = max(abs(rf.signal)/sys.gamma);     % Gauss
fwrite(fid, round(amp*C.RFSCALE), 'int16');

% waveform
rho = 2*round(abs(rf.signal/sys.gamma)/amp*C.MAXIAMP/2);  % int16
theta = 2*round(angle(rf.signal)/pi*C.MAXIAMP/2);         % int16
fwrite(fid, numel(rho), 'int16');   % number of samples in waveform
fwrite(fid, rho, 'int16');
fwrite(fid, theta, 'int16');

return


function sub_writegrad(fid, g, sys)

C = pulsegeq.constants;

% type
if isempty(g)
    fwrite(fid, C.NULL, 'int16');
    return;
end
if strcmp(g.type, 'trap')
    fwrite(fid, C.TRAP, 'int16');
else
    fwrite(fid, C.ARBITRARY, 'int16');
end

% delay
fwrite(fid, round(g.delay*1e6), 'int16');     % us

% amplitude
amp = g.amplitude/sys.gamma/100;  % Gauss/cm
fwrite(fid, round(amp*C.GSCALE), 'int16');

% waveform
if strcmp(g.type, 'trap')
    fwrite(fid, round(g.riseTime*1e6), 'int16');  % us
    fwrite(fid, round(g.flatTime*1e6), 'int16');
    fwrite(fid, round(g.fallTime*1e6), 'int16');
else
    % TODO: handle arbitrary gradients 
    % Remember to resample to 4us raster for ge
end

return


function sub_writeadc(fid, adc, sys)

C = pulsegeq.constants;

if isempty(adc)
    fwrite(fid, C.NULL, 'int16');
    return
end
fwrite(fid, C.ADC, 'int16');

fwrite(fid, adc.numSamples, 'int16');   
fwrite(fid, round(adc.dwell*1e6), 'int16');   % us
fwrite(fid, round(adc.delay*1e6), 'int16');   % us

return
