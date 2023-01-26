
% design sinc pulse
alpha = 30;       % degrees
slThick = 5e-3;   % m
[rf30_sinc, gz] = mr.makeSincPulse(alpha*pi/180, 'Duration',3e-3,...
    'SliceThickness', slThick, 'apodization', 0.42, 'timeBwProduct', 4, 'system',sys);

% simulate RF pulse (no gradients)
[M_z, M_xy, F2] = mr.simRf(rf30_sinc);

[bw, f0, M_xy_sta, F1] = mr.calcRfBandwidth(rf30_sinc);

figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, flip angle 30°');


% simulate RF pulse and gradients
figure;
raster = (rf30_sinc.t(2)-rf30_sinc.t(1));          % s
gamma = sys.gamma*1e-4;                            % Hz/Gauss
rf_g = rf.signal/gamma;                            % Gauss
rf_g = [zeros(1, rf.delay/raster) rf_g];
gz_gcm = pulsegeq.sub_trap2shape(gz, raster);         % Gauss/cm
rf_g = [rf_g zeros(1, length(gz_gcm)-length(rf_g))];  % make same length

m0 = [0 0 1];   % initial magnetization along mz
Z = 1.2*linspace(-slThick*100, slThick*100, 1000);         % cm
T1 = 1000; T2 = 100;                                       % ms
toppe.utils.rf.slicesim(m0, rf_g, gz_gcm, raster*1e3, Z, T1, T2);


