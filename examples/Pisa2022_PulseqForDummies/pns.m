close all

getsys;

sysGE.maxSlew = 20;  % G/cm/ms
sysGE.maxGrad = 5;  % G/cm

mxs = 20;
mxg = 5;

% create a gradient trapezoid spoiler
ncycles = 2;     % cycles of gradient spoiling across slthick
slthick = 0.5;   % cm
g1 = toppe.utils.makecrusher(ncycles,slthick,sysGE, 0, mxs,mxg);

g = g1;
toppe.writemod(sysGE, 'gz', g, 'ofname', 'tmp1.mod');
toppe.plotmod('tmp1.mod', 'gradient', 'xrm', 'printPNS', true);


% make it a multipolar gradient
g = [g1; -g1];
g = [g; g; g];
toppe.writemod(sysGE, 'gz', g, 'ofname', 'tmp3.mod');
toppe.plotmod('tmp3.mod', 'gradient', 'xrm', 'printPNS', false);

% create 'infinitely' long ramp to check stimulation threshold
sysGE.maxSlew = 20;  % G/cm/ms
sysGE.maxGrad = 50;  % G/cm
mxs = 7;  % close to stimulation threshold
mxg = sysGE.maxGrad;
dg = mxs*sysGE.raster*1e3;  % max gradient change per sample (per axis)
g = 0:dg:mxg;
g = [g fliplr(g) 0 0]';
toppe.writemod(sysGE, 'gz', g, 'ofname', 'tmp4.mod');
toppe.plotmod('tmp4.mod', 'gradient', 'xrm', 'printPNS', false);
