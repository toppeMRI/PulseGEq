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
(in addition to the 'official' Pulseq toolbox, `+mr`).

Working with Pulseq on GE scanners in the way described here requires the
[PulseGEq](https://github.com/toppeMRI/PulseGEq)
MATLAB toolbox, which in turn requires the
[TOPPE](https://github.com/toppeMRI/toppe)
MATLAB toolbox.


## Setup

### Get the source code (Linux)

To install the code needed to run the GE portion of this course, do:
```
$ git clone git@github.com:toppeMRI/PulseGEq.git
$ git clone git@github.com:toppeMRI/toppe.git
$ cd PulseGEq; git checkout develop
```

Add the `PulseGEq` and `toppe` directories to your MATLAB path.

The code in this demo can be found in:
```
$ cd PulseGEq/examples/Pisa2022_PulseqForDummies
```


## GE (TOPPE) scan file structure and scan instructions

gre.tar contains multiple files that work together to define the execution of the MR sequence.
A detailed description of these files can be found
[here](https://github.com/toppeMRI/toppe/blob/main/Files.md).
For a brief overview, see slide 2 of the 
[slide deck](https://docs.google.com/presentation/d/1YsY_6vehFSnzAYg_EjwvA8p9BDq8BgcDjDCcX4gc5BE/edit?usp=sharing)
for this presentation.

Note: to display the sequence in MATLAB, only the following files are needed:
```
modules.txt
scanloop.txt
module (.mod) files 
```

Detailed scanning instructions can be found 
[here](https://github.com/jfnielsen/TOPPEpsdSourceCode/) (private repository -- restricted to GE users).

The acquired data is saved in either Pfiles and/or ScanArchive files, as with any vendor or custom sequence.
The function `toppe.utils.loadpfile` may be used to load P-files.



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

Convert gre.seq to the 'TOPPE' file format:
```
>> getsys  % or:
>> sysGE = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
        'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
        'gradient', 'xrm');
>> pulsegeq.seq2ge('gre.seq', sysGE, 'verbose', true, 'tarFile', 'gre.seq.tar');
```

Display the GE sequence:
```
>> % The following functions read modules.txt, scanloop.txt, and .mod files from working folder
>> toppe.plotseq(1, 8, sysGE);  
>> nModsPerTR = 4;
>> toppe.playseq(nModsperTR, sysGE, 'nTRskip', nModsPerTR, 'gmax', 3, 'rhomax', 0.01);
```

**Current limitations:**
1. A gap (200-400us) is inserted after each Pulseq block.
2. Within each block, the gradient waveforms on each axis must start and end with zero.



## Example 2: Pulseq to GE conversion (2D GRE), more efficient implementation

First clean up TOPPE scan files from last example (remove .mod files, modules.txt, scanloop.txt).

Then create a version of the previous sequence that contains fewer blocks:
```
>> writeGradientEcho_4ge;             % creates gre_4ge.seq
```

Let's compare with the original sequence:
```
>> example2;   % or:
>> seq.read('gre.seq');
>> seq.plot('timeRange', [0 12e-3], 'showblocks', true);
>> seq.read('gre_4ge.seq');
>> seq.plot('timeRange', [0 12e-3], 'showblocks', true);
```

Convert gre_4ge.seq to the TOPPE file format:
```
>> pulsegeq.seq2ge('gre_4ge.seq', sysGE, 'verbose', true, 'tarFile', 'gre_4ge.seq.tar');
```

Display the GE sequence:
```
>> toppe.plotseq(1, 4, sysGE);
>> nModsPerTR = 2;
>> toppe.playseq(nModsperTR, sysGE, 'nTRskip', nModsPerTR, 'gmax', 3, 'rhomax', 0.01);
```


## Example 3: GE to Pulseq conversion (3D FLASH/SPGR/T1-FFE)

For sites that use GE as their primary development platform, 
and that wish to also export the sequence to Siemens scanners,
it may be a good option to first create a set of TOPPE files using the TOPPE MATLAB toolbox, 
and then convert the sequence to Pulseq using `ge2seq.m`.

For a detailed walk-through of pulse programming in TOPPE, see the 
[slide deck](https://docs.google.com/presentation/d/1YsY_6vehFSnzAYg_EjwvA8p9BDq8BgcDjDCcX4gc5BE/edit?usp=sharing)
for this presentation.


### Create the TOPPE sequence files (flash3d.tar)

First clean up TOPPE scan files from last example (remove .mod files, modules.txt, scanloop.txt).

Then:
```
>> writeflash3d_4ge;   % creates flash3d.tar
```

Display it:
```
>> toppe.plotseq(1, 4, sysGE, 'gmax', 5, 'rhomax', 0.04);
>> nModsPerTR = 2;
>> toppe.playseq(nModsperTR, sysGE, 'nTRskip', nModsPerTR, 'gmax', 3, 'rhomax', 0.04);
```

(Code walk-through)


### Convert flash3d.tar to flash3d.seq

```
>> example3;  % or:
>> pulsegeq.ge2seq('flash3d.tar', sys, sys, ...
    'seqFile', 'flash3d.seq', ...  % output file name
    'nt', 200, ...
    'FOV', FOV/100);                % m
>> seq = mr.Sequence(sys);
>> seq.read('flash3d.seq');
>> seq.plot('timeRange', [0 20e-3]);
```



