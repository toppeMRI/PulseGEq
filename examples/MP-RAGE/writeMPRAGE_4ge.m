% set system limits (slew rate 130 and max_grad 30 work on Prisma)
sys = mr.opts('MaxGrad', 24, 'GradUnit', 'mT/m', ...
    'MaxSlew', 100, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq = mr.Sequence(sys);           % Create a new sequence object
fov = [192 240 256]*1e-3;         % Define FOV and resolution
N = [192 240 256];                % matrix sizes
alpha = 7;                        % flip angle
%ro_dur=5017.6e-6; % BW=200Hz/pix
load GE   % see bw.m
ro_dur = GE.ro_dur;
%ro_os=2;                        % readout oversampling
ro_os = GE.decimation;
ro_spoil = 3;                    % additional k-max excursion for RO spoiling
TI = 1.1;
TRout=2.5;
% TE & TR in the inner loop are as short as possible derived from the above parameters and the system specs
% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rfLen=100e-6;
ax=struct; % encoding axes
ax.d1='z'; % the fastest dimension (readout)
ax.d2='y'; % the second-fast dimension (the inner pe loop)
ax.d3=setdiff('xyz',[ax.d1 ax.d2]); % automatically set the slowest dimension
ax.n1=strfind('xyz',ax.d1);
ax.n2=strfind('xyz',ax.d2);
ax.n3=strfind('xyz',ax.d3);

% Create alpha-degree hard pulse 
rf = mr.makeBlockPulse(alpha*pi/180, sys, 'Duration',rfLen);

% Create adiabatic inversion pulse
rf180 = mr.makeAdiabaticPulse('hypsec', sys, 'Duration', 10.24e-3, 'dwell',1e-5);

% Define other gradients and ADC events
deltak=1./fov;
gro = mr.makeTrapezoid(ax.d1, ...
    'Amplitude', N(ax.n1)*deltak(ax.n1)/ro_dur, ...
    'FlatTime', ceil(ro_dur/sys.gradRasterTime)*sys.gradRasterTime, ...
    'system',sys);
adc = mr.makeAdc(N(ax.n1)*ro_os, ...
    'Duration', ro_dur,...
    'Delay', gro.riseTime, ...
    'system',sys);
groPre = mr.makeTrapezoid(ax.d1, ...
    'Area', -gro.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gro.riseTime),...
    'system',sys); % the first 0.5 is necessary to acount for the Siemens sampling in the center of the dwell periods
gpe1 = mr.makeTrapezoid(ax.d2, ...   
    'Area', deltak(ax.n2)*(N(ax.n2)/2),...  
    'system',sys);  % maximum PE1 gradient
gpe2 = mr.makeTrapezoid(ax.d3, ...
    'Area', deltak(ax.n3)*(N(ax.n3)/2), ...
    'system',sys);  % maximum PE2 gradient
gslSp = mr.makeTrapezoid(ax.d3, ...
    'Area', max(deltak.*N)*4, ...  % spoil with 4x cycles per voxel
    'Duration', 10e-3, ...
    'system',sys);  

% we cut the RO gradient into two parts for the optimal spoiler timing
[gro1,groSp]=mr.splitGradientAt(gro,gro.riseTime+gro.flatTime);
% gradient spoiling
if ro_spoil>0
    groSp = mr.makeExtendedTrapezoidArea(gro.channel, gro.amplitude, 0, deltak(ax.n1)/2*N(ax.n1)*ro_spoil, sys);
end

% calculate timing of the fast loop 
% we will have two blocks in the inner loop:
% 1: RF 
% 2: prewinder,phase enconding + readout + spoilers/rewinders
[groPre] = mr.align('right',groPre,mr.makeDelay(mr.calcDuration(gpe1,gpe2)-gro.riseTime));
gro1.delay = mr.calcDuration(groPre);
groSp.delay = mr.calcDuration(gro1);
adc.delay = gro1.delay+gro.riseTime;
gro1 = mr.addGradients({gro1,groPre,groSp},'system',sys);
gpe1c = mr.addGradients({gpe1, ...
    mr.makeTrapezoid(ax.d2, 'Area', -gpe1.area, 'duration', groSp.shape_dur, 'delay', groSp.delay, 'system',sys)});
gpe2c = mr.addGradients({gpe2, ...
    mr.makeTrapezoid(ax.d3, 'Area', -gpe2.area, 'duration', groSp.shape_dur, 'delay', groSp.delay, 'system',sys)});
TRinner = mr.calcDuration(rf)+mr.calcDuration(gro1); % we'll need it for the TI delay

% peSteps -- control reordering
pe1Steps = ((0:N(ax.n2)-1)-N(ax.n2)/2)/N(ax.n2)*2;
pe2Steps = ((0:N(ax.n3)-1)-N(ax.n3)/2)/N(ax.n3)*2;

% TI calc
TIdelay=round((TI-(find(pe1Steps==0)-1)*TRinner-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay)-rf.delay-mr.calcRfCenter(rf))/sys.blockDurationRaster)*sys.blockDurationRaster;
TRoutDelay=TRout-TRinner*N(ax.n2)-TIdelay-mr.calcDuration(rf180);

% pre-register objects that do not change while looping
gslSp.id=seq.registerGradEvent(gslSp);
gro1.id=seq.registerGradEvent(gro1);
[~, gpe1c.shapeIDs]=seq.registerGradEvent(gpe1c);
[~, gpe2c.shapeIDs]=seq.registerGradEvent(gpe2c);
[~, rf.shapeIDs]=seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes 
[rf180.id, rf180.shapeIDs]=seq.registerRfEvent(rf180); % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create .mod files and modules.txt needed for execution on GE scanners

gam = 4.2576e3;     % Hz/Gauss

% set system specs
GE.sys = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...  % <= scanner limit
    'maxGrad', 7, 'gradUnit', 'Gauss/cm', ...  % maxGrad *MUST* match the physical scanner spec 
    'maxRf', 0.25, 'rfUnit', 'Gauss'); % <= scanner limit

% set TOPPE low-level system timing variables (scanner specific; times are in us)
% overrides defaults
GE.sys.toppe.version = 'v3';
GE.sys.toppe.start_core_daq = 100;
GE.sys.toppe.myrfdel = 148;
GE.sys.toppe.daqdel = 156;
GE.sys.toppe.timetrwait = 64;
GE.sys.timessi = 200;

% write inversion pulse to inversion.mod
GE.rf180.signal = [0; rf180.signal(:)/gam; 0];  % Gauss
GE.rf180.signal = toppe.utils.makeGElength(GE.rf180.signal);  % enforce 4-sample boundary
toppe.writemod('rf', GE.rf180.signal, 'ofname', 'inversion.mod');

% write spoil.mod
GE.spoil.mxs = 8;  % Gauss/cm/ms. Lower to reduce PNS.
res = fov(1)/N(1)*100;   % spatial resolution (cm)
GE.spoil.nSpoilCycles = 4;
gCrush = toppe.utils.makecrusher(GE.spoil.nSpoilCycles, res, 0, GE.spoil.mxs, GE.sys.maxGrad);
toppe.writemod('gz', gCrush, 'ofname', 'spoil.mod', 'system', GE.sys);

% write alpha pulse to tipdown.mod
GE.rf.n = 10;         % number of RF samples (will be padded with zeros below)
GE.rf.dur = GE.rf.n*GE.raster; 
GE.rf.signal = [0; (alpha/360) / (gam * GE.rf.dur) * ones(GE.rf.n,1); 0]; % must start and end with zeros
GE.rf.signal = toppe.utils.makeGElength(GE.rf.signal);  % enforce 4-sample boundary
toppe.writemod('rf', GE.rf.signal, 'ofname', 'tipdown.mod');

% create readout.mod for TOPPE
GE.zres = fov(ax.n3)/N(ax.n3) * 1e2;   % cm
[GE.ro.wav, GE.pe1.wav, GE.pe2.wav] = toppe.utils.makegre(fov(ax.n1)*1e2, N(ax.n1), GE.zres, ...
    'oprbw', GE.oprbw, ...
    'ncycles', 2, ...    % number of cycles of spoiling along readout (x)
    'system', GE.sys, ...
    'ofname', 'tmp.mod');
toppe.writemod('gx', GE.pe2.wav, 'gy', GE.pe1.wav, 'gz', GE.ro.wav, ...
    'system', GE.sys, 'ofname', 'readout.mod');
system('rm tmp.mod');

% create modules.txt file for TOPPE
% If placing in a folder other than /usr/g/bin/,
% you must provide the full filename path.
GE.modFileText = ['' ...
'Total number of unique cores\n' ...
'4\n' ...
'fname  duration(us)    hasRF?  hasDAQ?\n' ...
'inversion.mod\t0\t1\t0\n' ...   % entries are tab-separated   
'spoil.mod\t0\t0\t0\n' ...
'tipdown.mod\t0\t1\t0\n' ...
'readout.mod\t0\t0\t1' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, GE.modFileText);
fclose(fid);

%% done creating .mod files and modules.txt for TOPPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define GRAPPA undersampling pattern along slowest dimension (outer pe loop)
nACS = 20;  % number of auto-calibration lines (dense sampling in center)
R = 2;      % undersampling outside ACS region
S = 0*(1:N(ax.n3));
S(1:R:end) = 1;
%S( (end/2-nACS/2):(end/2+nACS/2-1) ) = 1;  % TODO: uncomment
J = find(S == 1);

% intialize scanloop.txt file for TOPPE
toppe.write2loop('setup', 'version', 3);

% Scan loop
for i = 1:N(ax.n2)  % add noise scans for the .seq file
    seq.addBlock(adc);  
end
GE.slice = 2;   % reserve slice = 1 for noise scans (at end of scan)
for j = J  % 1:N(ax.n3) 
    % inversion pulse, spoiler, and delay
    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(TIdelay),gslSp);

    % for TOPPE
  	toppe.write2loop('inversion.mod', ...
		'RFamplitude', 1.0);
	toppe.write2loop('spoil.mod', ...
		'textra', round(TIdelay*1e3)); % (ms) approximate TODO

    rf_phase=0;
    rf_inc=0;

    % pre-register the PE gradients that repeat in the inner loop
    gpe2cj=mr.scaleGrad(gpe2c,pe2Steps(j));
    gpe2cj.id=seq.registerGradEvent(gpe2cj);

    for i=1:N(ax.n2)
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);

        % excitation and readout (for .seq file)
        seq.addBlock(rf);
        seq.addBlock(adc,gro1,mr.scaleGrad(gpe1c,pe1Steps(i)),gpe2cj);

        % excitation and readout (for TOPPE)
        toppe.write2loop('tipdown.mod', ...
            'RFphase', rf_phase/180*pi); 
        GE.textra = (i == N(ax.n2)) * round(TRoutDelay*1e3); % approximate TODO
        toppe.write2loop('readout.mod', ...
            'Gamplitude', [pe2Steps(j) pe1Steps(i) 1.0]', ...
            'DAQphase', rf_phase/180*pi, ...
            'textra', GE.textra, ...
            'slice', GE.slice, 'view', i);
    end
    
    GE.slice = GE.slice + 1;

    seq.addBlock(mr.makeDelay(TRoutDelay));
end

% add noise scans for TOPPE
% do this last since receive gain is based on signal from beginning of sequence
% first insert ~5s pause to allow magnetization/system to settle
toppe.write2loop('spoil.mod', ...
    'Gamplitude', [0 0 0]', ...
    'textra', 5e3);
for i = 1:N(ax.n2)  
    toppe.write2loop('readout.mod', ...
        'Gamplitude', [0 0 0]', ... % turn off gradients 
        'slice', 1, 'view', i);   % GE data is stored in 'slice', 'echo', and 'view' indeces
end

% finalize scanloop.txt
toppe.write2loop('finish');   

fprintf('Sequence ready\n');

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% plot, etc
seq.plot('TimeRange',[0 TRout*2]);

%%
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'mprage');

seq.write('mprage.seq')       % Write to pulseq file
%seq.install('siemens');

return;

%% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
%if Nx<=32
    tic;
    [kfa,ta,kf]=seq.calculateKspacePP();
    toc
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
%end



%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 



