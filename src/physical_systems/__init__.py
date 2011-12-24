"""
The physical system package contains classes that represent
different quantum mechanical systems.

A physical system should at least define a B{potential} and the
degrees of freedom (B{DOF}) in the system. The potential should be a string
and be given in a format that is understood by OpenCL. This is
very similar to C however some "native" functions can be used such
as native_exp.

The variables in the potential are given as xn where n is an integer
representing the particular degree of freedom. An harmonic oscillator
in three dimensions with spring constant 1.0 would thus be given as: 

M{0.5f * (x1 * x1 + x2 * x2 + x3 * x3)}

Eventual operators that are defined uses the same way of writing. 
"""
