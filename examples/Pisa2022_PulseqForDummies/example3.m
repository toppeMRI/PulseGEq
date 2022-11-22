
pulsegeq.ge2seq('flash3d.tar', sysGE, sys, ...
    'seqFile', 'flash3d.seq', ...  % output file name
    'nt', 100, ...                 % only convert the first 100 blocks (for testing)
    'FOV', FOV/100);               % m

seq = mr.Sequence(sys);
seq.read('flash3d.seq');
seq.plot('timeRange', [0 20e-3]);

