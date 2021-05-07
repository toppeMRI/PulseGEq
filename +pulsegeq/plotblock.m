function plotblock(blk, gradRasterTime, rfRasterTime)
% function plotblock(blk, gradRasterTime, rfRasterTime)
%
% Inputs:
%  blk     Pulseq block, e.g., blk = seq.getBlock(10);
%  grad/rf raster time:  sec

gamma = 4257.6; % Hz/G

gx = [];
gy = [];
gz = [];
rf = [];

tgx = 0;
tgy = 0;
tgz = 0;
trf =0;

subplot(121); hold on; 
xlabel('ms');
ylabel('G/cm');

if ~isempty(blk.gx)
	gx = pulsegeq.sub_trap2shape(blk.gx, gamma, gradRasterTime);
	tgx = linspace(blk.gx.delay+gradRasterTime, blk.gx.delay+gradRasterTime*length(gx), length(gx));
	plot(tgx, gx, 'r.'); 
end
if ~isempty(blk.gy)
	gy = pulsegeq.sub_trap2shape(blk.gy, gamma, gradRasterTime);
	tgy = linspace(blk.gy.delay+gradRasterTime, blk.gy.delay+gradRasterTime*length(gy), length(gy));
	plot(tgy, gy, 'g.'); 
end
if ~isempty(blk.gz)
	gz = pulsegeq.sub_trap2shape(blk.gz, gamma, gradRasterTime);
	tgz = linspace(blk.gz.delay+gradRasterTime, blk.gz.delay+gradRasterTime*length(gz), length(gz));
	plot(tgz, gz, 'b.'); 
end
if ~isempty(blk.rf)
   %rf = downsample(blk.rf.signal/gamma, round(raster/1e-6));  % Gauss
	rf = blk.rf.signal/gamma; % Gauss
	trf  = linspace(blk.rf.delay+rfRasterTime, blk.rf.delay+rfRasterTime*length(rf), length(rf));
	subplot(122); 
	plot(trf, rf, 'b.'); ylabel('rf'); xlabel('ms');
	legend('rf (G)');
end

tbeg = min([tgx(1) tgy(1) tgz(1) trf(1)]);
tend = max([tgx(end) tgy(end) tgz(end) trf(end)]);

subplot(121); 
legend('gx (G/cm)', 'gy', 'gz');
gmax = max([max(gx) max(gy) max(gz)]);
gmin = min([min(gx) min(gy) min(gz)]);
if ~isempty(gx) | ~isempty(gy) | ~isempty(gz)
	axis([max(0,tbeg-0.5e-3) tend+0.5e-3 gmin-0.1 gmax+0.1]);
end

subplot(122);
if ~isempty(rf)
	axis([max(0,tbeg-0.5e-3) tend+0.5e-3 min(rf)-0.05 max(rf)+0.05]);
end


