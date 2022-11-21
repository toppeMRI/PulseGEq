# Learning with examples: Part II

"Pulseq course for dummies" satellite course, Pisa, Italy 23-Nov-2022


## Sequence analysis: slice profile

With Pulseq, we have complete knowledge of the MR sequence, which
allows us to inspect and check the sequence before exporting it to the scanner.

Here we will simulate the slice profile, to
verify that the flip angle and slice thickness are correct.
Knowledge of the (theoretical) slice profile is also useful in quantitative MRI settings.

sinc pulse:
```
>> rfsim;
```

Simultaneous multi-slice (SMS) pulse:
```
>> rfsim_sms;
```


## Sequence analysis: peripheral nerve stimulation (PNS)

GE model + implementation


## Sequence analysis: specific absorption rate (SAR)




## Non-cartesian imaging example: 3D stack of spirals

stack of spirals, including fat saturation

simulate spectral profile of fat sat pulse
