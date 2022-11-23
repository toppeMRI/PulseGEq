# Learning with examples: Part II

"Pulseq course for dummies" satellite course, Pisa, Italy 23-Nov-2022

With Pulseq, we have complete knowledge of the MR sequence, which
allows us to **inspect and check the sequence** before exporting it to the scanner.



## Sequence analysis: peripheral nerve stimulation (PNS)

GE's model is described in Schulte et al 2014.

Example:
```
>> pns;
```

A related collection of papers can be found [here](https://www.sciencedirect.com/topics/neuroscience/chronaxie).

Information about Siemens' SAFE model can be found
[here](https://github.com/filip-szczepankiewicz/safe_pns_prediction).



## Sequence analysis: RF pulses


### slice profile

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


### RF power (SAR) and gradient heating

SAR can be estimated in advance (offline) by comparing the RF power 
to a reference sequence with known SAR monitor reading for a certain patient weight.

Similarly, the gradient power can be calculated offline and entered into the
interpreter's gradient heating checks.

GE example:
`toppe.preflightcheck()`

See also the [sar4seq](https://github.com/imr-framework/sar4seq) Github repository (Geethanath et al).

