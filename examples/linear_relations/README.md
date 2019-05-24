Linear relations
================

In this example, we show how to detect linear relations between the
output elements of a dataflow graph, before performing their analytic
reconstruction.  In particular, we use this in order to write the most
complex elements in terms of simpler ones, thanks to the possibility
of estimating the complexity of multivariate graphs via simpler
univariate reconstruction.  We then proceed to reconstruct a
complete independent subset of output functions, which we use to give
a more compact form of the result.

The example contains the following files:

- [full_output.m](full_output.m) contains a representation of the full
  output of the graph, which is not meant to be read or used directly.
- [linear_relations.wl](linear_relations.wl) computes the linear
  relations and the functional reconstruction discussed above.
