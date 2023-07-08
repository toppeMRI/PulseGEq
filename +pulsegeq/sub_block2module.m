function module = sub_block2module(block, nextblock, blockid, system, modnum)
%% Storage for groups of arbitrary gradients and rf signals
%
% Initialize a module with waveforms from 'block'
%
% Inputs:
%   block            Pulseq block obtained with getBlock()
%   nextblock        Pulseq block obtained with getBlock(). Needed in case it's a delay block.
%   blockid          Pulseq block is seq.getBlock(blockid)
%   system           hardware specs, as described in ../seq2ge.m
%   modnum           Module number. Determines .mod file name.
% Output:
%   module      Struct containing module info (rf/gradient waveforms, blockEvents, etc)
%               This corresponds to the TP_MODULE struct in the TOPPE interpreter
%                 .duration        'duration' entry in modules.txt
%                 .hasRF           'hasRF' entry in modules.txt
%                 .hasADC          'hasADC' entry in modules.txt
%                 .fname           .mod output file name
%                 .blockids        
%                 .modnum          .mod file id/index (order in modules.txt)
%                 .rf              [nt npulses] Normalized RF waveforms 
%                 .gx              [nt npulses] Normalized Gx waveforms (same for .gy, .gz)
%                 .gy 
%                 .gz 
%                 .npulses         number of different sets of waveforms (e.g., size(gx.waveforms,2))
%                 .npres  
%                 .rfres  
%                 .res  

import pulsegeq.*

%if length(moduleArr)+1 > system.toppe.nMaxModules
%	error(sprintf('The number of modules exceeds system.toppe.nMaxModules (%d).\nAre you sure you need that many modules?', system.toppe.nMaxModules));
%end

% Initialize with defaults
module          = struct();
module.duration = block.blockDuration*1e6;
module.hasRF    = 0;
module.hasADC   = 0;
module.fname   = sprintf('module%d.mod', modnum);
module.blockids = [];
module.modnum   = modnum;
	
% Initialize waveform arrays to []
for ax = {'rf', 'gx','gy','gz'}
	module.(ax{1}) = [];
end

module.npulses = 0;

% Update 'module' with waveform information from 'block'.
module = sub_updatemodule(module, block, blockid, system);  

return
