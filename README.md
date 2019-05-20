FiniteFlow MathTools
====================

This is a collection of packages, tools and examples using the
Mathematica interface to the FiniteFlow library.

The examples presented here are designed to be simple enough to run in
a few minutes on a modern laptop.  They also contain plenty of
comments which should be useful as an introduction to using FiniteFlow
for solving these types of problems.  However, they can be used as
templates for more complicated applications to be run on larger
machines and clusters.

It should be noted that the packages included in this repository
should be reguarded as a set of utilities rather the implementation of
fully automated solutions for specific tasks.

The algorithms implemented here, are based on the following
publication

- Tiziano Peraro, *FiniteFlow: multivariate functional reconstruction
  using finite fields and dataflow graphs* (to appear)


Installation and usage
----------------------

In order to use and run these packages and examples, follow the steps
below:

* Install the [FiniteFlow](https://github.com/peraro/finiteflow)
  library, following the instructions in its`README.md` file
* The LiteIBP package in this repository and some of the examples
depend on [LiteRed](http://www.inp.nsk.su/~lee/programs/LiteRed/).
You can download the latter from the previous link.  In order to make
sure its files are in your Mathematica path, consider adding the
following to your `init.m` file
```
$LiteRedPath = "/path/to/litered"
If[Not[MemberQ[$Path,$LiteRedPath]],$Path = Flatten[{$Path, $LiteRedPath }]];
```
You can skip this step if you are not interested in the LiteIBP
package or in the examples using it.
* Make sure the packages in this repository are in your Mathamatica
path, e.g. by adding the following to your `init.m` file
```
$FFMathTools = "/path/to/this/repo/packages"
If[Not[MemberQ[$Path,$FFMathTools]],$Path = Flatten[{$Path, $FFMathTools }]];
```

Most of the examples must be run from their own directory.


Contents
--------

This reposiotory contains the following subdirectories:

* [packages/](packages/README.md): contains a collection of Mathematica
  packages and tools which use the FiniteFlow library
* [examples/](examples/README.md): contains a collection of examples which use
  the FiniteFlow library and the Mathematica packages in this
  repository


Known issues
------------

Some Mathematica functions in the LiteIBP package may not work in
Mathematica 11.3 and later if any parallel kernel is launched in the
main program (see Mathematica's `Kernels[]` and `LaunchKernels[]`).
The issue only involves Mathematica scripts generating and serializing
IBPs, and it does not affect parallelization in the functional
reconstruction implemented in FiniteFlow.  In order to prevent this
issue, we recommend using the "LaunchKernels" option of the LiteIBP
procedure, as in the examples of this repostory.
