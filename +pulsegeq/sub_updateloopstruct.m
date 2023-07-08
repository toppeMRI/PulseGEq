%% Update/initialize loopStructArr 
function arg = sub_updateloopstruct(arg, block, system, varargin)
%
% Fill struct containing entries for one row in scanloop.txt 
%
% Inputs
%  arg          Either [] (empty), or a loopStruct struct (see below)
%  block        Pulseq block (getBlock()). Can be empty.
%  system       See ../seq2ge.m
%
% Options
%  see code below

import pulsegeq.*

% Initialize loopStruct struct
if isempty(arg)
    % Defaults
    arg.mod   = 1;           % module number. Positive integer (starts at 1).
    arg.rfamp = 0;           % rf waveform (rho) amplitude (Gauss)
    arg.gxamp = 0;           % Gx waveform amplitude       (Gauss/cm)
    arg.gyamp = 0;           % Gy waveform amplitude
    arg.gzamp = 0;           % Gz waveform amplitude
    arg.slice = 1;           % data storage ’slice’ index. Positive integer (starts at 1).
    arg.echo  = 0;           % data storage ’echo’ index.  Non-negative integer (starts at 0).
    arg.view  = 1;           % data storage ’view’ index.  Positive integer (starts at 1).
    arg.dabmode = 1;         % turn on/off data acquisition (1/0)
    arg.phi    = 0;          % in-plane (x-y) rotation angle (radians)
    arg.rfphs  = 0;          % RF transmit phase (radians)
    arg.recphs = 0;          % receive phase (radians)
    arg.textra = 0;          % time added to end of module (sec)
    arg.rffreq = 0;          % RF transmit frequency offset (Hz)
    arg.wavnum = 1;          % waveform number (rf/grad waveform array column index). Non-zero positive integer.
    arg.rotmat = eye(3);     % 3x3 rotation matrix (added to toppev3)
    arg.blockGroupID = [];
end

% Substitute specified system values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

% If block is provided, fill in values as appropriate
if ~isempty(block)
    % Set block group id
    if isfield(block, 'label')
        arg.blockGroupID = block.label.value;
    else
        arg.blockGroupID = 0;
    end

    % If pure delay block
    if pulsegeq.isdelayblock(block)
        minWaitPulseDuration = 4e-6;   % s
        arg.mod = 0;  % tells the interpreter to play wait pulse
        arg.textra = block.blockDuration - minWaitPulseDuration;
    end

    % rf amplitude, phase, freq
    if ~isempty(block.rf)
        arg.rfamp = max(abs(block.rf.signal)/system.gamma);    % Gauss
        arg.rfphs = block.rf.phaseOffset;                      % radians
        arg.rffreq = round(block.rf.freqOffset);                      % Hz
    end

    % receive phase
    if ~isempty(block.adc) 
        arg.recphs = block.adc.phaseOffset;                      % radians
    end

    % gradient amplitude
    for ax = {'gx','gy','gz'};
        ax = cell2mat(ax);
        grad = block.(ax);
        if ~isempty(grad)
            if strcmp(grad.type,'trap')
                % Gauss/cm, signed (in modules array, traps are normalized and positive)
                eval( sprintf('arg.%samp = (block.%s.amplitude)/system.gamma/100;', ax, ax) );   
            else
                % arbitrary waveform
                eval( sprintf('wav = block.%s.waveform/system.gamma/100;', ax) );    % Gauss/cm
                wmax = max(abs(wav));
                eval( sprintf('arg.%samp = wmax;', ax) );
            end
        else
        end
    end
end

return;

