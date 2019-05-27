# BioEnergetics

This repository contains the combined bioenergetics and visual
foraging model used for the [GrowChinook
website](http://growchinook.fw.oregonstate.edu/). It is
packaged as a python module which may be called independently of the
GrowChinook website.

## Installation

To install this Python package, clone this repository to your
computer, navigate to its root directory (where `setup.py` is
located), and type:

```pip install .```

Alternatively, you can install it directly from github without cloning

```pip install git+git://github.com/reservoirwebs/bioenergetics```

### Dependencies

This library has the following dependencies. Either of the above
installation methods will install them automatically.

- Numpy
- Scipy
- Matplotlib

## Documentation

API documentation in HTML format can be found in the [doc
directory](./doc), and is hosted
[here](http://growchinook.fw.oregonstate.edu/apidocs/).

Have a look in the [unit tests](./tests) directory for example usage.
