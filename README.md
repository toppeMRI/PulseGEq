# PulseGEq

Pulseq (http://pulseq.github.io/) is an open, vendor-agnostic file format specification for MR pulse sequences.
Pulseq files can be created in a number of ways, e.g., using the Pulseq Matlab package `+mr`, or 
[JEMRIS](http://jemris.org/) (additional [notes](JEMRIS.md)).

The code in this repository converts a Pulseq (.seq) file to a set of files that can be executed on GE scanners.

The Pulseq and TOPPE repositories are included here as Git submodules (./deps/),
to enure that the correct Pulseq/TOPPE versions (commits) are used.


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
  >> addpath ~/github/toppeMRI/PulseGEq/deps/pulseq/matlab/    % +mr package
  >> addpath ~/github/toppeMRI/PulseGEq/deps/toppe/            % +toppe package
```


## Usage example

```
  >> cd ./examples/
  >> pulsegeq.seq2ge('2DFLASH.seq', 'verbose', true);
```
or
```
  >> seq = mr.Sequence();
  >> seq.read('2DFLASH.seq');
  >> seq.plot('timeRange', [0 0.04]);
  >> pulsegeq.seq2ge(seq, 'verbose', true);
```

To display sequence:
```
  >> nModsPerTR = 3;
  >> toppe.playseq(nModsPerTR);

```

Screen capture of this example: https://www.youtube.com/embed/qswI1vPQ4io


## Support for older Pulseq versions


### Pulseq v1.2.1

As of Oct 5, 2020, JEMRIS outputs Pulseq v1.2.1.
To use that version of Pulseq, do:


#### Set your clone of the Pulseq repo to v1.2.1

```
$ git clone git@github.com:pulseq/pulseq.git
$ git checkout 74eb4c06d66ca60a6a6d8548d3ccc1584bca0b98
```
The clone is now in a detached state. Later, to get back to the branch you were on, do
```
$ git checkout -
```

#### Call seq2ge.m as follows
```
pulsegeq.seq2ge('gre.seq', 'pulseqVersion', 'v1.2.1');
```
