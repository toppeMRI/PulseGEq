% Convert 2DFLASH.seq to a set of files that can be executed
% on GE scanners.
%
% For an explanation of these files, see
% https://github.com/toppeMRI/toppe/blob/main/Files.md
%
% For scan instructions using these files, see
% https://github.com/jfnielsen/TOPPEpsdSourceCode/

% Create a subdirectory where we'll work
system('mkdir -p sandbox');
cd('sandbox');

% get the various Matlab packages and add to path
system('git clone git@github.com:pulseq/pulseq.git');
system('git clone git@github.com:toppeMRI/toppe.git');
system('git clone git@github.com:toppeMRI/PulseGEq.git');

addpath('pulseq/matlab');   % +mr package
addpath('toppe');           % +toppe package
addpath('PulseGEq');        % +pulsegeq package

% Create *.mod files, scanloop.txt, and modules.txt from 2DFLASH.seq
sys = toppe.systemspecs('maxSlew', 20, 'maxGrad', 50, 'gradient', 'xrm');
pulsegeq.seq2ge('../2DFLASH.seq', sys, 'verbose', false);
system('rm toppeScanFiles.tar');

% For TOPPE v4 we also need a 'seqstamp.txt' file
toppe.preflightcheck('../toppe0.entry', 'seqstamp.txt', sys);

% Place the file ../toppe0.entry in /usr/g/research/pulseq/ on the scanner,
% and all other files in /usr/g/research/pulseq/2dflash/.
% Then scan as described in https://github.com/jfnielsen/TOPPEpsdSourceCode/README.md

% Optional steps
system('tar cf scanfiles.tar ../toppe0.entry *.mod modules.txt scanloop.txt seqstamp.txt');

toppe.plotseq(1, 10, sys);
toppe.plotmod('all', 'gradcoil', sys.gradient);
nModsPerTR = 3;
figure; toppe.playseq(nModsPerTR, sys);



