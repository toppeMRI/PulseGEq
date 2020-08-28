# PulseGEq

Pulseq (http://pulseq.github.io/) is an open, vendor-agnostic file format specification for MR pulse sequences.

The code in this repository converts a Pulseq (.seq) file to a set of files that can be executed on GE scanners.

The Pulseq repository is included here as a Git submodule (lib/pulseq).


## Get the source code

To ensure that the submodule(s) are loaded, add the `--recursive` option to the git clone command:

<!--- 
$ git clone --recurse-submodules git@github.com:toppeMRI/PulseGEq.git>
-->

```
  $ git clone --recursive git@github.com:toppeMRI/PulseGEq.git
```


## Set Matlab paths

```
  >> addpath ~/github/toppeMRI/PulseGEq/src/                  % seg2ge.m
  >> addpath ~/github/toppeMRI/PulseGEq/lib/pulseq/matlab/    % +mr package
  >> addpath ~/github/toppeMRI/PulseGEq/lib/toppe/            % +toppe package
```


## Usage example

```
  >> cd ../examples/
  >> seq2ge('2DFLASH.seq', 'verbose', true);
```
or
```
  >> seq = mr.Sequence();
  >> seq.read('2DFLASH.seq');
  >> seq.plot('timeRange', [0 0.04]);
  >> seq2ge(seq, 'verbose', true);
```

To display sequence:
```
  >> nModsPerTR = 3;
  >> toppe.playseq(nModsPerTR);

```

Screen capture of this example: https://www.youtube.com/embed/qswI1vPQ4io
