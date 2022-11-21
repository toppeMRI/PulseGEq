# Learning with examples: Part II

"Pulseq course for dummies" satellite course, Pisa, Italy 23-Nov-2022


## Sequence analysis: slice profile

With Pulseq, we have complete knowledge of the MR sequence, which
allows us to inspect and check the sequence before exporting it to the scanner.

Here we will simulate the slice profile, to
verify that the flip angle and slice thickness are correct.
Knowledge of the (theoretical) slice profile is also useful in quantitative MRI settings.

```
alpha = 30;
[rf, gz] = mr.makeSincPulse(alpha*pi/180, 'Duration',3e-3,...
    'SliceThickness', 5e-3, 'apodization', 0.42, 'timeBwProduct', 4, 'system',sys);

[M_z, M_xy, F2] = mr.simRf(rf30_sinc);

[bw, f0, M_xy_sta, F1] = mr.calcRfBandwidth(rf30_sinc);
```


### Design the slice-selective excitation

```
>> rf30_sinc = mr.makeSincPulse(pi/6,'system',sys,'Duration',3e-3,'use','excitation',...
    'PhaseOffset',pi/2,'apodization',0.3,'timeBwProduct', 8);
>> [M_z,M_xy,F2] = mr.simRf(rf30_sinc);
>> [bw,f0,M_xy_sta,F1] = mr.calcRfBandwidth(rf30_sinc);
```


+mr function

toppe function

sinc pulse example
SLR pulse example

SMS pulse example (w/ imaging results)


## Sequence analysis: peripheral nerve stimulation (PNS)

Vendors measure PNS experimentally to determine threshold.

GE model + implementation


## Sequence analysis: specific absorption rate (SAR)




## Non-cartesian imaging example: 3D stack of spirals

stack of spirals, including fat saturation

simulate spectral profile of fat sat pulse
