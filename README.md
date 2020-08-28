# PulseGEq

Pulseq (http://pulseq.github.io/) is an open, vendor-agnostic file format specification for MR pulse sequences.

The code in this repository converts a Pulseq (.seq) file to a set of files that can be executed on GE scanners.

The Pulseq repository is included here as a Git submodule (lib/pulseq).


## Get the source code

To ensure that the submodule(s) are loaded, add the `--recurse-submodules` option to the git clone command:

<!--- 
$ git clone --recurse-submodules git@github.com:toppeMRI/PulseGEq.git>
-->

```
  $ git clone --recursive git@github.com:toppeMRI/PulseGEq.git
```

