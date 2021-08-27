% readout bandwidth calculation

% must match writeMPRAGE_4ge.m
fov = 25.6;  % cm
nx = 256;

% 'oprbw' control variable (CV); bandwidth/fov is 2*oprbw (kHz)
oprbw = 125/5;

% design readout gradient
tmp = 1;  % z resolution (cm) (dummy)
ofname = 'readout.mod';
toppe.utils.makegre(fov, nx, tmp, ...
    'oprbw', oprbw, ...
    'ofname', ofname);
[rf,gx,gy,gz,desc,hdr] = toppe.readmod(ofname);
gx = gx((hdr(1)+1):(hdr(1)+hdr(2)));  % get readout plateau portion

% BW per pixel
bwpp = 2*oprbw/nx*1e3;   % Hz

% save readout parameters
GE.raster = 4e-6;   % (s) ADC dwell time, gradient/RF raster, are all the same
GE.decimation = 125/oprbw;  % readout oversampling factor
GE.ro_dur = GE.raster * nx * GE.decimation;
GE.fov = fov;
GE.nx = nx;
GE.oprbw = oprbw;

fprintf('GE readout parameters:\n   raster = %.3e s; ro_dur = %.3e; nsamples = %d;\n   grad amp = %.6f G/cm, bw/pix = %.3f Hz\n', ...
    GE.raster, GE.ro_dur, nx * GE.decimation, gx(1), bwpp);

save GE GE

%ro_dur=5017.6e-6; % BW=200Hz/pix

