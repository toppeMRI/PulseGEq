function writecores(fid, cores)
% cores = cell array, with each entry containing a variable-length vector
%         containing the integer IDs of the blocks in one core:
%         [<core ID> <nBlocksInCore> <block ID 1> <block ID 2> ...]

nCores = length(cores);

fwrite(fid, nCores, 'int16');

for ic = 1:nCores
    fwrite(fid, cores{ic}(1), 'int16');  % core ID
    fwrite(fid, cores{ic}(2), 'int16');  % number of blocks in core
    fwrite(fid, cores{ic}(3:end), 'int16');  % block IDs
end
    
