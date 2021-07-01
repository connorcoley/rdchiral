[![PyPI version](https://badge.fury.io/py/rdchiral.svg)](https://badge.fury.io/py/rdchiral)

# rdchiral
Wrapper for RDKit's RunReactants to improve stereochemistry handling

## Requirements

* RDKit (version >= 2019)
* Python (version >= 3.5)

## Installation

To install RDChiral run

```pip install rdchiral```

To get the most recent version reflected by this git repo, install with

```pip install -e "git://github.com/connorcoley/rdchiral.git#egg=rdchiral"```

## Fast C++ version (rdchiral_cpp)

A fast C++ implementation ([rdchiral_cpp](https://gitlab.com/ljn917/rdchiral_cpp)) is available as a drop-in replacement. It provides ~10x speedup. To install from anaconda, run

```conda install -c conda-forge -c ljn917 rdchiral_cpp```

## Documentation

See ```rdchiral/main.py``` for a brief description of expected behavior and a few basic examples of how to use the wrapper. 

See ```rdchiral/test/test_rdchiral.py``` for a small set of test cases described [here](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00286)
