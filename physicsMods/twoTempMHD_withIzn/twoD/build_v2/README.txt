this is the third version with full 2D set of equations.
The new additions are electron heat flux and ion viscosity.
The heat flux has magnetic field corrections. The viscosity 
is unmagnetized. The maximum normalized ion collision time
used for the viscosity coefficient is set in the input deck.
I should probably do something similar for electron heat flux.

I'm using relaxation scheme for both heat flux and viscosity
since these can be diffusive in nature. I currently only
have the xx component, which is the only non-zero component in 1D
done using relaxtion, I need to do the other two as well.

Maybe define maximum mean free paths for heat flux and viscosity?
