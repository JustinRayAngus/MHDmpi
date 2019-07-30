build_v0:
This is initial version with full 2D set of equations.
It does not include finite-ion or finite-electron 
inertial length corrections. 

build_v1:
This version included finite-ion and finite-electron
inertial length corrections, though it doesn't seeem
to be working all that well. Can't decide between using
J or J0 when doing all these extra corrections. I think
I need to do proper upwinding to pressure term in Ohm's
law. Try turning ionization off. Look at why I have so
much noise-like stuff in front of shock front. What parts
of E/J are these from?

build_v2:
This version has electron heat flux and ion viscosity.
Using relaxation scheme for both.
