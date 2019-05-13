Differential equations
======================

An example computing differential equations for a basis of UT master
integrals.  The topology is a double box with one off-shell leg, but
the example can be easily adapted to more complex topologies.

The example contains the following files:

- [topology.wl](topology.wl) contains the definition of the topology
  and the external kinematics.  It is imported by the other files and
  it should not be run directly.
- [generate_ibps.wl](generate_ibps.wl]) generates the system of IBP
  identities and serializes it in JSON format.  The definition of the
  UT integrals is also added to the system as additional equations.
  This file must be run from this directory before the ones below.
- [differential_equations.wl](differential_equations.wl) reconstructs
  the differential equations for the master integrals, without making
  any assumption on their functional form.  It must be run from this
  directory after `generate_ibps.wl`.
- [differential_equations_from_ansatz.wl](differential_equations_from_ansatz.wl)
  reconstructs the differential equations for the master integrals,
  from an ansatz based on the alphabet of the topology, which is
  assumed to be known apriori.  It must be run from this directory
  after `generate_ibps.wl`.
