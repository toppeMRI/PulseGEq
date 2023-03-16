function Cores = readcores(fid)
% Cores = cell array, with each entry containing a variable-length vector
%         containing the integer IDs of the blocks in one core:
%         [<core ID> <nBlocksInCore> <block ID 1> <block ID 2> ...]

nCores = fread(fid, 1, 'int16');

for ic = 1:nCores
    coreID = fread(fid, 1, 'int16');
    nBlocksInCore = fread(fid, 1, 'int16');
    blockIDs = fread(fid, nBlocksInCore, 'int16');
    Cores{ic} = [coreID nBlocksInCore blockIDs'];
end
    
