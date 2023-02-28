function writecores(fid, cores)
% cores = cell array, with each entry containing a variable-length vector
%         [<core ID> <nBlocksInCore> <core ID 1> <core ID 2> ...]

nCores = length(cores);

fwrite(fid, nCores, 'int16');

for ic = 1:nCores
    fwrite(fid, cores{ic}, 'int16');
end
    
