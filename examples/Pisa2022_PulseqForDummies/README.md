# Pulseq for dummies (GE related content)

## Installing the source code

At present, working with Pulseq on GE scanners requires the
[PulseGEq](https://github.com/toppeMRI/PulseGEq)
MATLAB toolbox, which in turn requires the
[TOPPE](https://github.com/toppeMRI/toppe)
MATLAB toolbox.

To install the code needed to run the GE portion of this course, do:
```
$ mkdir temp
$ cd temp
$ git clone  TODO (Pulseq, ...)
```

## Enter the demo directory and set MATLAB path

```
$ cd PulseGEq/examples/Pisa2022_PulseqForDummies
```

```
>> setup;
```

## Creating a 2D FLASH sequence 

### Create the .seq file (2dflash.seq)
```
>> write2DFLASH;
```

### Check and plot .seq file


## Convert 2dflash.seq file to the 'TOPPE' file format

### Set hardware limits 

```
sys.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ...     % psd_rf_wait (gradient/rf delay, us)
    'daqdel', 152, ...      % psd_grd_wait (gradient/acquisition delay, us)
    'gradient', 'hrmb');    % xrm: MR750; hrmb: UHP; hrmw: Premier
```

### Do the file conversion
```
>>  pulsegeq.seq2ge('2dflash.seq', sys.ge, 'verbose', false);
```
