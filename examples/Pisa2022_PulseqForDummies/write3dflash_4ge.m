function write3dflash_4ge

sysGE = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ...                          % psd_rf_wait (gradient/rf delay, us)
    'daqdel', 152, ...                           % psd_grd_wait (gradient/acquisition delay, us)
    'gradient', 'xrm');                          % xrm: MR750; hrmb: UHP; hrmw: Premier

N = [100 100 50];   % matrix size
FOV = [20 20 10];   % FOV (cm)
flip = 10;          % degrees
DTE = [0];          % extend TE by this much (ms). If vector, acquisitions are interleaved

if false
sub_flash(sysGE, N, FOV, flip, DTE, ...
    'rfSpoilSeed', 117, ...
    'nCyclesSpoil', 2, ...
    'fatsat', false);
end

rf = toppe.readmod('tipdown.mod');
toppe.plotseq(1, 2, sysGE, ...
    'rhomax', 1.1*max(abs(rf))); 

return

function sub_flash(sys, N, FOV, flip, DTE, varargin)
% function sub_flash(sys, N, FOV, flip, DTE, varargin)
%
% Fully-sampled 3D RF-spoiled GRE sequence for B0 (and B1-) mapping.
% 
% Inputs:
%  sys        struct   scanner hardware settings. See toppe.systemspecs().
%  N          [1 3]    matrix size
%  FOV        [1 3]    field of view (cm)
%  flip       [1 1]    flip angle (degrees)
%  DTE        [1 n]    increase (from minimum value ) in TE for each of n scans (ms)
%
% Input options with defaults:
%  entryFile = 'toppeN.entry';
%  scanFilePath = '/usr/g/research/pulseq/cal/b0/';  
%  tbw = 8;                     % RF pulse time-bandwidth product
%  rfDur = 2;                   % RF pulse duration (ms)
%  ftype = 'min';               % 'min': minimum-phase SLR pulse; 'ls': linear phase
%  slabThick = 0.8*FOV(3);      % excited slab thickness
%  rfSpoilSeed = 117;              % RF spoiling phase increment factor (degrees)
%  exMod        = 'tipdown.mod';
%  readoutMod   = 'readout.mod';
%  nCyclesSpoil = 2;               % number of cycles of phase across voxel (along x and z)
%  fatsat       = false;           % struct containing fat sat settings:
%  fatFreqSign = -1;               %

% defaults
arg.entryFile = 'toppeN.entry';
arg.scanFilePath = '/usr/g/research/pulseq/cal/b0/';
arg.tbw = 8;                     % RF pulse time-bandwidth product
arg.rfDur = 2;                   % RF pulse duration (ms)
arg.ftype = 'min';               % 'min': minimum-phase SLR pulse; 'ls': linear phase
arg.slabThick = 0.8*FOV(3);      % excited slab thickness
arg.rfSpoilSeed = 117;           % RF spoiling phase increment factor (degrees)
arg.exMod         = 'tipdown.mod';
arg.readoutMod    = 'readout.mod';
arg.nCyclesSpoil = 2;   % number of cycles of phase across voxel (along x and z)
arg.fatsat       = false;         % add fat saturation pulse?
arg.fatFreqSign = -1;            % sign of fatsat pulse frequency offset

% substitute with provided keyword arguments
arg = toppe.utils.vararg_pair(arg, varargin);

nScans = numel(DTE);

% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if N(1) ~= N(2) | FOV(1) ~= FOV(2)
    error('In-plane FOV and matrix be square.');
end

voxSize = FOV./N;  % cm

ny = N(2);
nz = N(3);

% Write modules.txt
nModules = 2 + arg.fatsat;
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', nModules);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
if arg.fatsat
    fprintf(fid, '%s\t0\t1\t0\n', 'fatsat.mod');
end
fprintf(fid, '%s\t0\t1\t0\n', arg.exMod);
fprintf(fid, '%s\t0\t0\t1\n', arg.readoutMod);
fclose(fid);

% Write entry file.
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.scanFilePath, ...
    'b1ScalingFile', arg.exMod, ...
    'readoutFile', arg.readoutMod);


%% Create .mod files

if arg.fatsat
    % fat sat module
    fatsat.flip    = 90;
    fatsat.slThick = 1000;       % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
    fatsat.tbw     = 2.0;        % time-bandwidth product
    fatsat.dur     = 4.5;        % pulse duration (ms)

    b1 = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, sys, ...
        'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
        'writeModFile', false);
    b1 = toppe.makeGElength(b1);
    toppe.writemod(sys, 'rf', b1, 'ofname', 'fatsat.mod', 'desc', 'fat sat pulse');
end

% excitation module
[ex.rf, ex.g] = toppe.utils.rf.makeslr(flip, arg.slabThick, ...
    arg.tbw, arg.rfDur, nz*arg.nCyclesSpoil, sys, ...
    'ftype', arg.ftype, ...
    'spoilDerate', 0.5, ...
    'ofname', arg.exMod);

% readout module
% Here we use the helper function 'makegre' to do that, but that's not a requirement.
% Reduce slew to keep PNS in normal mode (<80% of limit)
toppe.utils.makegre(FOV(1), N(1), voxSize(3), sys, ... 
    'ofname', arg.readoutMod, ...
    'ncycles', arg.nCyclesSpoil); 


%% Write scanloop.txt
rfphs = 0;              % radians
rfSpoilSeed_cnt = 0;
ny = N(2);
nz = N(3);

toppe.write2loop('setup', sys, 'version', 4);  % initialize file ('scanloop.txt')

for iz = -1:nz     % We use iz<1 for approach to steady-state
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%d of %d', max(1,iz), nz);
    for iy = 1:ny
        for ite = 1:length(DTE)
            % y/z phase encode amplitudes. Turn off during approach to steady-state.
            % My convention is to start at (-kymax, -kzmax)
            a_gy = -((iy-1+0.5)-ny/2)/(ny/2) * (iz>0);  
            a_gz = -((iz-1+0.5)-nz/2)/(nz/2) * (iz>0);

            if(arg.fatsat)
                fatChemShift = 3.5;  % fat/water chemical shift (ppm)
                fatFreq = arg.fatFreqSign*sys.gamma*1e4*sys.B0*fatChemShift*1e-6;  % Hz
                toppe.write2loop('fatsat.mod', sys, ...
                    'RFoffset', round(fatFreq), ...   % Hz
                    'RFphase', rfphs);         % radians

                rfphs = rfphs + (arg.rfSpoilSeed/180*pi)*rfSpoilSeed_cnt ;  % radians
                rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
            end

            toppe.write2loop(arg.exMod, sys, ...
                'RFamplitude', 1.0, ...
                'textra', DTE(ite), ...
                'RFphase', rfphs);

            toppe.write2loop(arg.readoutMod, sys, ...
                'Gamplitude', [1.0 a_gy a_gz]', ...
                'DAQphase', rfphs, ...
                'textra', max(DTE) - DTE(ite), ... % to keep TR constant
                'slice', max(iz,1), 'echo', ite, 'view', iy);

            rfphs = rfphs + (arg.rfSpoilSeed/180*pi)*rfSpoilSeed_cnt ;  % radians
            rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);  % finalize file

fprintf('TR = %.3f ms\n', toppe.getTRtime(1, 2, sys)*1e3);

%figure; toppe.plotseq(1, 4, sys);

% Create 'sequence stamp' file for TOPPE
% This file is listed in line 6 of the .entry file
toppe.preflightcheck(arg.entryFile, 'seqstamp.txt', sys);

% Write files to tar archive (for convenience).
system(sprintf('tar cf b0.tar %s seqstamp.txt scanloop.txt modules.txt *.mod', arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, sys, ...
%    'tpause', 0.05, ...
%    'nTRskip', 8);

return;

