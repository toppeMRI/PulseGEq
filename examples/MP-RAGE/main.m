addpath ~/github/pulseq/matlab/       % path to Pulseq +mr Matlab package
addpath ~/github/toppeMRI/PulseGEq/   % path to +pulsegeq Matlab package
addpath ~/gitlab/fMRI/toppe/          % path to +toppe Matlab package

% create .seq file. Blocks are chosen to make the seq2ge.m conversion easy
writeMPRAGE_4ge;  


