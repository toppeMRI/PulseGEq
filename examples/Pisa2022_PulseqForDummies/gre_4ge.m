
sysGE = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ...                          % psd_rf_wait (gradient/rf delay, us)
    'daqdel', 152, ...                           % psd_grd_wait (gradient/acquisition delay, us)
    'gradient', 'xrm');                          % xrm: MR750; hrmb: UHP; hrmw: Premier

pulsegeq.seq2ge('gre.seq', sysGE, 'verbose', true, 'nt', 100);
