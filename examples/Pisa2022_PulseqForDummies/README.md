# Pulseq for dummies (GE related content)

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

### Create the .seq file (2dgre.seq)
```
>> write2dgre;
```

==> Check and plot .seq file


### Convert 2dgre.seq file to the 'TOPPE' file format

Set GE scanner hardware limits 
```
sys.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ...                          % psd_rf_wait (gradient/rf delay, us)
    'daqdel', 152, ...                           % psd_grd_wait (gradient/acquisition delay, us)
    'gradient', 'xrm');                          % xrm: MR750; hrmb: UHP; hrmw: Premier
```

Do the file conversion
```
>>  pulsegeq.seq2ge('2dgre.seq', sys.ge, 'verbose', false);
```


### Execute the sequence on a GE scanner


### Load the data and reconstruct


## Example 2: GE to Pulseq conversion (3D FLASH/SPGR/T1-FFE)

For sites that use GE as their primary development platform, 
and that wish to also export the sequence to Siemens scanners,
it may be a good option to first create a set of TOPPE files using the TOPPE MATLAB toolbox, 
and then convert the sequence to Pulseq using `ge2seq.m`.


### Create the TOPPE sequence files (3dflash.tar)

```
>> write3dflash_4ge;
```

The file 3dflash.tar can be executed on GE scanners as described above for the 2D GRE scan.

3dflash.tar contains multiple files that work together to define the execution of the MR sequence.
A detailed description of these files can be found
[here](https://github.com/toppeMRI/toppe/blob/main/Files.md).


### Convert 3dflash.tar to 3dflash.seq

Set Siemens scanner hardware limits
```
>> sys.siemens = mr.
```

Convert to .seq file
```
>> pulsegeq.ge2seq('3dflash.tar', sys.siemens, 'ofname', '3dflash.seq');
```


