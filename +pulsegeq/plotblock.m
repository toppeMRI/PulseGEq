function plotblock(blk, gamma, raster)
% 3T: gamma = 4257.6 Hz/G
% raster: 4e-6; 

rf = [];
gx = [];
gy = [];
gz = [];

if ~isempty(blk.gx)
	gx = pulsegeq.sub_trap2shape(blk.gx, gamma, raster);
end
if ~isempty(blk.gy)
	gy = pulsegeq.sub_trap2shape(blk.gy, gamma, raster);
end
if ~isempty(blk.gz)
	gz = pulsegeq.sub_trap2shape(blk.gz, gamma, raster);
end
if ~isempty(blk.rf)
   rf = downsample(blk.rf.signal/gamma, round(raster/1e-6));  % Gauss
end
	
subplot(121); hold on; 
plot(raster*1e3*(1:length(gx)), gx, 'r.'); 
plot(raster*1e3*(1:length(gy)), gy, 'g.'); 
plot(raster*1e3*(1:length(gz)), gz, 'b.'); 
xlabel('ms');
ylabel('G/cm');
legend('gx (G/cm)', 'gy', 'gz');
if ~isempty(rf)
	subplot(122); 
	plot(raster*1e3*(1:length(rf)), rf); ylabel('rf'); xlabel('ms');
	legend('rf (G)');
end

