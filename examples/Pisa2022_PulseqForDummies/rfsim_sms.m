% load pulse that excites multiple slices
%

% load RF pulse and gradients
[rf_g, ~, ~, gz_gcm] = toppe.readmod('tipdown_sms.mod.git');

% simulate RF pulse and gradients
figure;
raster = 4e-6;   % sec
m0 = [0 0 1];    % initial magnetization along mz (a.u.)
fov = 8;                                      % FOV for simulation (cm)
Z = 1.2*linspace(-fov/2, fov/2, 1000);        % cm
T1 = 1000; T2 = 100;                                       % ms
toppe.utils.rf.slicesim(m0, rf_g, gz_gcm, raster*1e3, Z, T1, T2);


