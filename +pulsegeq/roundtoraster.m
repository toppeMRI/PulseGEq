function dur = roundtoraster(dur, raster)

dur = dur + raster - mod(dur, raster);

