function test(seqFile)

seq = mr.Sequence();
seq.read(seqFile);

sys = pulsegeq.systemspecs('maxSlew', 20, 'maxGrad', 5);

% test writeblock.m and readblock.m
if false
    b = seq.getBlock(1);
    fid = fopen('block.bin', 'w', 'ieee-be');
    pulsegeq.writeblock(fid, 7, b, sys);
    fclose(fid);

    fid = fopen('block.bin', 'r', 'ieee-be');
    b2 = pulsegeq.readblock(fid);
    fclose(fid);

    b
    b2
end

% test writeGEseq.m and readGEseq.m
if true
    pulsegeq.seq2ge(seqFile, sys, 'nt', 20);  % creates out.4ge
    [ParentBlocks, Cores, Dyn] = pulsegeq.readGEseq('out.4ge');
    b = seq.getBlock(1);
    ParentBlocks{1};
    keyboard
end


