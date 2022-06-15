# Conductances

A `ConductanceSystem` describes the state- and time-dependence of a conductance (i.e.
classically _g_, the inverse of resistance; typically measured in Siemens). By default, a
`Conudctance System` is composed of zero or more gates and a parameter for the maximum
conductance value (``\overline{g}``). The product of each gate's output and ``\overline{g}``
yield absolute conductance over time.

`ConductanceSystem` can be used to describe the conductance associated with ionic membrane currents, synaptic currents, and axial currents that flow between connected neuronal compartments. 

```math
g(t) = \overline{g}m^3h
```

!!! note
    We intentionally separate the dynamics  current and conductance. 
