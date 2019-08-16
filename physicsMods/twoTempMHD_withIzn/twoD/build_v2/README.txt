this is the third version with full 2D set of equations.
The new additions are electron heat flux and ion viscosity.
The heat flux has magnetic field corrections. The viscosity 
is unmagnetized. The maximum normalized ion collision time
used for the viscosity coefficient is set in the input deck.
I should probably do something similar for electron heat flux.

I'm using relaxation scheme for both heat flux and viscosity
since these can be diffusive in nature.

Maybe define maximum mean free paths for heat flux and viscosity?


Aug 16, 2019

Signicicant refactoring and additions. Added full stress tensor
to Ohm's law and electron viscosity. The viscosity is working, but
only set for 1D right now. Need to revisit the tensor and divergence
calculations considering velocity at cell edges. Added both of these
to try and mitigate small scale structures in J formed when using
finite ion scale corrections. I have not tried ion-scale corrections
since getting viscosity to work because I need to get viscosity
working for full 2D first.


