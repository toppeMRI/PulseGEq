% Create readout.mod and tipdown.mod, to replace module4.mod and module3.mod,
% respectively (as created by pulsegeq.seq2ge('mprage.seq').
% This is done to try to match timing.

sys = toppe.systemspecs;

% sequence parameters
load GE

% create readout.mod
ofname = 'readout.mod';
tmp = 1;  % z resolution (cm) (dummy)
toppe.utils.makegre(GE.fov, GE.nx, tmp, ...
    'oprbw', oprbw, ...
    'ncycles', 2, ...    % number of cycles of spoiling along readout (x)
    'system', sys, ...
    'ofname', ofname);

