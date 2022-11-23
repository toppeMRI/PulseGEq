close all

getsys;

sysGE.maxSlew = 20;  % G/cm/ms
sysGE.maxGrad = 5;  % G/cm

mxs = 15;
mxg = 5;

% create a gradient trapezoid spoiler
ncycles = 2;     % cycles of gradient spoiling across slthick
slthick = 0.5;   % cm
g = toppe.utils.makecrusher(ncycles,slthick,sysGE, 0, mxs,mxg);

toppe.writemod(sysGE, 'gx', g, 'gz', g, 'ofname', 'tmp1.mod');
toppe.plotmod('tmp1.mod', 'gradient', 'xrm', 'printPNS', true);

return

% compare with a bipolar gradient
g = [g; -g];
toppe.writemod(sysGE, 'gx', g, 'gz', g, 'ofname', 'tmp2.mod');
toppe.plotmod('tmp2.mod', 'gradient', 'xrm', 'printPNS', false);

% create 'infinitely' long ramp
dg = mxs
toppe.plotmod('tmp2.mod', 'gradient', 'test', 'printPNS', false);
