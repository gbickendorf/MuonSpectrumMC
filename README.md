# MuonSpectrumMC

This project is an expansion of my master thesis. The goal is to  derive 
the additional contribution to the muon decay spectrum due to an 
additional field that serves as a dark matter candidate. The 
corresponding matrix elements were calculated in mathematica and 
exported in C++ syntax.
This time the integration was carried out by a Monte Carlo integrator 
found in the GNU Scientific Library. 
To extract bounds on the parameter space, this spectrum was superimposed
with the expected standard model spectrum. Using the simpex minimisation
technique of GSL on the squared-difference function to the expected 
behavior could then be used to extract the experimentally well known 
Michel-Parameters. By finding the coupling constant that saturates the
experimental limits we could extract constraints on the parameter space.


Changes from the MuonSpectrumCUDA project:
- Chance quadrature integration on the GPU to GSL Monte Carlo integration
- Add analysis tools that were before left for other software.
