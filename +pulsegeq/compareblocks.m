function issame = compareblocks(b1, b2)
    % Compare two Pulseq blocks

    issame = true;

    if b1.blockDuration ~= b2.blockDuration
        issame = false; return;
    end
    if (xor(isempty(b1.rf), isempty(b2.rf)) | ...
        xor(isempty(b1.gx), isempty(b2.gx)) | ... 
        xor(isempty(b1.gy), isempty(b2.gy)) | ... 
        xor(isempty(b1.gz), isempty(b2.gz)) | ... 
        xor(isempty(b1.adc), isempty(b2.adc)))
        issame = false; return;
    end
    if ~isempty(b1.rf)
        if ~comparerf(b1.rf, b2.rf)
            issame = false; return;
        end
    end
    if ~isempty(b1.gx)
        if ~comparegradients(b1.gx, b2.gx)
            issame = false; return;
        end
    end
    if ~isempty(b1.gy)
        if ~comparegradients(b1.gy, b2.gy)
            issame = false; return;
        end
    end
    if ~isempty(b1.gz)
        if ~comparegradients(b1.gz, b2.gz)
            issame = false; return;
        end
    end
    if ~isempty(b1.adc)
        if ~compareadc(b1.adc, b2.adc)
            issame = false; return;
        end
    end
return

function issame = comparerf(rf1, rf2)

    issame = true;

    if (rf1.delay ~= rf1.delay | ...
        rf1.shape_dur ~= rf2.shape_dur | ...
        norm(rf1.t-rf2.t) > 1e-6)
        issame = false; return;
    end

    % compare (normalized) signal shapes
    wav1 = abs(rf1.signal)/max(abs(rf1.signal));
    wav2 = abs(rf2.signal)/max(abs(rf2.signal));
    if norm(wav1 - wav2) > 1e-4
        issame = false; return;
    end

return

function issame = comparegradients(g1, g2)

    issame = true;

    if ~strcmp(g1.type, g2.type)
        issame = false; return;
    end
    if strcmp(g1.type, 'trap')
        if (g1.riseTime ~= g2.riseTime | ...
            g1.flatTime ~= g2.flatTime | ...
            g1.fallTime ~= g2.fallTime | ...
            g1.delay ~= g2.delay)
            issame = false; return;
        end
    else
        % TODO: handle arbitrary waveforms
        tol = 1e-4;
    end

return

function issame = compareadc(adc1, adc2)

    issame = true;

    if (adc1.numSamples ~= adc2.numSamples | ...
        adc1.dwell ~= adc2.dwell | ...
        adc1.delay ~= adc2.delay)
        issame = false;
    end

return

