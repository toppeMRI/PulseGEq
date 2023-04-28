function [ParentBlocks, Cores, Dyn] = readGEseq(fname)

C = pulsegeq.constants();

if C.LITTLEENDIAN
    fid = fopen(fname, 'r', 'ieee-le');
else
    fid = fopen(fname, 'r', 'ieee-be');
end

nParentBlocks = fread(fid, 1, 'int16')
for blockID = 1:nParentBlocks
    ParentBlocks{blockID} = pulsegeq.readblock(fid);
end
Cores = pulsegeq.readcores(fid);
Dyn = [];

% dynamics (scan loop)
n1 = fread(fid, 1, 'int16');
n2 = fread(fid, 1, 'int16');
nt = n1*C.MAXIAMP + n2

for ib = 1:nt
    Dyn(ib, :) = fread(fid, 9, 'int16');
end

fclose(fid);
