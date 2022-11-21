# Pulseq on GE

"Pulseq course for dummies" satellite course, Pisa, Italy 23-Nov-2022

## Overview

At present, Pulseq compatibility with GE relies on a **file conversion** step, 
whereby the .seq file is loaded into MATLAB and used to write a set of files
that can be executed on GE scanners using the TOPPE interpreter.
This conversion is not exact, but is sufficient for many sequence types.

The reverse is also possible: a set of TOPPE scan files can be used to generate the
corresponding .seq file.
That conversion is exact.
Thus, the TOPPE MATLAB toolbox can be viewed as an alternative set of tools for creating a .seq file
(in addition to the 'official' Pulseq toolbox).

Working with Pulseq on GE scanners requires the
[PulseGEq](https://github.com/toppeMRI/PulseGEq)
MATLAB toolbox, which in turn requires the
[TOPPE](https://github.com/toppeMRI/toppe)
MATLAB toolbox.


## Setup

### Get the source code

To install the code needed to run the GE portion of this course, do:
```
$ mkdir temp
$ cd temp
$ git clone  TODO (Pulseq, ...)
```

### Enter the demo directory and set MATLAB path

```
$ cd PulseGEq/examples/Pisa2022_PulseqForDummies
```

```
>> setup;    % add the +mr, +pulsegeq, and +toppe toolboxes to the MATLAB path
```


## Example 1: Pulseq to GE conversion (2D GRE)

Create the .seq file (gre.seq):
```
>> writeGradientEcho; 
```

Plot the .seq file:
```
>> seq = mr.Sequence(sys);
>> seq.read('gre.seq');
>> seq.plot('timeRange', [0 12e-3], 'showblocks', true);
```

Convert gre.seq file to the 'TOPPE' file format:
```
>> sysGE = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ...                          % psd_rf_wait (gradient/rf delay, us)
    'daqdel', 152, ...                           % psd_grd_wait (gradient/acquisition delay, us)
    'gradient', 'xrm');                          % xrm: MR750; hrmb: UHP; hrmw: Premier

>> pulsegeq.seq2ge('gre.seq', sysGE, 'verbose', true, 'tarFile', 'gre.seq.tar');
```

Display the GE sequence:
```
>> nModsPerTR = X;
>> toppe.playseq(nModsperTR, sysGE, 'nTRskip', nModsPerTR, 'gmax', 3, 'rhomax', 0.01);
```

## GE file structure

gre.tar contains multiple files that work together to define the execution of the MR sequence.
A detailed description of these files can be found
[here](https://github.com/toppeMRI/toppe/blob/main/Files.md).



## Example 2: Pulseq to GE conversion (2D GRE), more efficient implementation

Create a version of the previous sequence that contains fewer blocks:
```
>> system('rm *.mod');       % remove the GE 'module' files
>> writeGradientEcho_4ge;    % creates gre_4ge.seq
```

Convert gre_4ge.seq file to the TOPPE file format:
```
>> pulsegeq.seq2ge('gre_4gre.seq', sysGE, 'verbose', true, 'tarFile', 'gre_4ge.seq.tar');
```

Display the GE sequence:
```
>> nModsPerTR = 2;
>> toppe.playseq(nModsperTR, sysGE, 'nTRskip', nModsPerTR, 'gmax', 3, 'rhomax', 0.01);
```


## Executing the sequences on a GE scanner and loading the data




## Example 3: GE to Pulseq conversion (3D FLASH/SPGR/T1-FFE)

For sites that use GE as their primary development platform, 
and that wish to also export the sequence to Siemens scanners,
it may be a good option to first create a set of TOPPE files using the TOPPE MATLAB toolbox, 
and then convert the sequence to Pulseq using `ge2seq.m`.


### Create the TOPPE sequence files (flash3d.tar)

```
>> writeflash3d;
```

The file flash3d.tar can be executed on GE scanners as described above for the 2D GRE scan.

flash3d.tar contains multiple files that work together to define the execution of the MR sequence.
A detailed description of these files can be found
[here](https://github.com/toppeMRI/toppe/blob/main/Files.md).


### Convert flash3d.tar to flash3d.seq

(flash3d2seq.m)

```
sysSiemens = mr.opts('MaxGrad', 50, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

pulsegeq.ge2seq('flash3d.tar', sys, sysSiemens, ...
    'seqFile', 'flash3d.seq', ...  % output file name
    'nt', 200, ...
    'FOV', FOV/100);                % m
```



