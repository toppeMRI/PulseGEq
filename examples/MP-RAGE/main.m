addpath ~/github/pulseq/matlab/       % path to Pulseq +mr Matlab package
addpath ~/github/toppeMRI/PulseGEq/   % path to +pulsegeq Matlab package
addpath ~/gitlab/fMRI/toppe/          % path to +toppe Matlab package

% create GE.mat
bw;

% create .seq file. Blocks are chosen to make the seq2ge.m conversion easy
% Loads readout gradient (gx) parameters from GE.mat
writeMPRAGE_4ge;  


