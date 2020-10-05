# JEMRIS notes

JEMRIS website: http://jemris.org/


## Prerequisites

Ubuntu packages:
```
libsundials-cvode2   libsundials-dev
libxerces-c3.2  libxerces-c-dev
libhdf5-dev   libhdf5-cpp-100
libcln6   libcln-dev
libginac6   libginac-dev
```

“Note that all libraries above are mandatory for compiling JEMRIS (runtime and development versions need to be present)”
Compiling JEMRIS (pulseq-export branch)
Download the JEMRIS source package from Github. Change into the pulseq-export branch.

## Installation

Get source code
```
$ git clone git@github.com:JEMRIS/jemris.git
```

Compile and install
```
$ cd jemris
$ mkdir build; cd build; cmake ..
$ make
$ ctest -V
$ sudo make install
```

## Run the Sequence GUI and export a sequence to Pulseq

In Matlab:
```
>> addpath /usr/local/share/jemris/matlab
>> JEMRIS_seq;
```
Load the share/examples/gre.xml example sequence and ‘export to scanner’ to write a Pulseq (.seq) file


