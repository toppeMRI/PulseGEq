% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=256; Ny=256;     % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=3e-3;            % slice
TR=12e-3;                       % repetition time TR
TE=5e-3;                        % echo time TE  
%TE=[7.38 9.84]*1e-3;            % alternatively give a vector here to have multiple TEs (e.g. for field mapping)

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
% B0=2.89; % 1.5 2.89 3.0
% sat_ppm=-3.45;
% sat_freq=sat_ppm*1e-6*B0*lims.gamma;
% rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
%     'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
% gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',3.2e-3,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',1e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTE>=0));
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

rf_phase=0;
rf_inc=0;

% create 'composite' block containing pe encoding gradients, readout gradient, and spoiler gradient

% First create dummy sequence, then write waveforms as aribtrary shapes
c = 1;
seq_d = mr.Sequence(sys);           % dummy sequence
gyPre = mr.makeTrapezoid('y','Area', -phaseAreas(1), 'Duration', mr.calcDuration(gxPre), 'system',sys);
seq_d.addBlock(gxPre,gyPre,gzReph);
seq_d.addBlock(mr.makeDelay(delayTE(c)));
seq_d.addBlock(gx,adc);
gyPre.amplitude=-gyPre.amplitude;
seq_d.addBlock(mr.makeDelay(delayTR(c)), gxSpoil, gyPre, gzSpoil);

% make new ADC
adcDur=seq_d.duration-2*sys.adcDeadTime; % dead times at the beginning and at the end 
adcSamplesPerSegment=100; % we need some "roundish" number of samples
adcNumSam=floor(adcDur/adc.dwell/adcSamplesPerSegment)*adcSamplesPerSegment; 
adc_read = mr.makeAdc(adcNumSam,'dwell',adc.dwell,'system',sys);

wave_data = seq_d.waveforms_and_times();
wave_length = seq_d.duration/sys.gradRasterTime;
wave_time = ((1:wave_length)-0.5)*sys.gradRasterTime;
wave_x = interp1(wave_data{1}(1,:), wave_data{1}(2,:), wave_time,'linear',0); % the last value of 0 is important to extrapolate with 0s
wave_y = interp1(wave_data{2}(1,:),wave_data{2}(2,:),wave_time,'linear',0);
wave_z = interp1(wave_data{3}(1,:),wave_data{3}(2,:),wave_time,'linear',0);

gx_read = mr.makeArbitraryGrad('x', [0 wave_x 0]);
gy_read = mr.makeArbitraryGrad('y', [0 wave_y 0]);
gz_read = mr.makeArbitraryGrad('z', [0 wave_z 0]);

%seq.addBlock(rf,gz);
%seq.addBlock(gx_read, gy_read, gz_read, adc_read);

%rf = mr.makeArbitraryRf(ex.signal, alpha/180*pi, 'system', sys);

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    gy_read_amplitude_scale = ((i-1) - Ny/2)/(Ny/2);
    for c=1:length(TE)
        %seq.addBlock(rf_fs,gz_fs); % fat-sat
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        %
        seq.addBlock(rf,gz);
        gy_read = mr.makeArbitraryGrad('y', gy_read_amplitude_scale*wave_y);
        %gy_read.amplitude = gy_read.amplitude * gy_read_amplitude_scale;
        seq.addBlock(gx_read, gy_read, gz_read, adc_read);
        %gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
        %seq.addBlock(gxPre,gyPre,gzReph);
        %seq.addBlock(mr.makeDelay(delayTE(c)));
        %seq.addBlock(gx,adc);
        %gyPre.amplitude=-gyPre.amplitude;
        %seq.addBlock(mr.makeDelay(delayTR(c)),gxSpoil,gyPre,gzSpoil)
    end
end

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre');

seq.write('gre_4ge.seq')       % Write to pulseq file

%seq.install('siemens');

%% plot sequence

seq.plot('timeRange', [0 2]*TR);

return

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

