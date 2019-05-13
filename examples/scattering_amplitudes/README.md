Scattering amplitudes
=====================

An example of the reduction of a scattering amplitude to
polylogarithms.  We show how to perform the reconstruction of the
coefficients of the Laurent expansion of the final result.  In
particular, we do not reconstruct the analytic expression of
intermediate stages, such as IBP reduction tables.

This is a simple one-loop four-point amplitude, but the example can be
adapted to more complex cases.  In particular, we follow a strategy
which is also suitable for multi-loop amplitudes and form factors.

The example contains the following files:

- [topology.wl](topology.wl) contains the definition of the topology
  and the external kinematics.  It is imported by the other files and
  it should not be run directly.
- [generate_ibps.wl](generate_ibps.wl]) generates a system of IBP
  identities for the amplitude and serializes it in JSON format.  This
  file must be run from this directory before the ones below.
- [unreduced_amplitude_-+-+.m](unreduced_amplitude_-+-+.m) is an
  unreduced helicity amplitude in terms of standard Feynman integrals.
  It is used by the example
  [compute_amplitude_simpler.wl](compute_amplitude_simpler.wl) and it
  is not meant to be run directly.
- [numerator.m](numerator.m) is the numerator of the amplitude in
  terms of loop momenta and polarization vectors, obtained using
  Feynman diagrams and Feynman rules.  It is used by the example
  [compute_amplitude.wl](compute_amplitude.wl) and it is not meant to
  be run directly.
- [compute_amplitude_simpler.wl](compute_amplitude_simpler.wl)
  reconstructs an helicity amplitude in terms of polylogs, starting
  from a representation of it in terms of (unreduced) standard Feynman
  integrals.  We recommend reading this example before the following
  one.
- [compute_amplitude.wl](compute_amplitude.wl) reconstructs an
  helicity amplitude in terms of polylogs, starting from a
  representation of it in terms of loop momenta and polarization
  vectors, obtained by substituting the Feynman rules into the
  diagrams.  This example uses integrand reduction and transverse
  integration before the IBP reduction, and it is thus more
  complicated than the previous one, but more useful when adapted to
  more complex calculations.

