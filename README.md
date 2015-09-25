# pykida

Python interface for the KInetic Database for Astrochemistry (KIDA) and the associated Nahoon code.

One of the issues we ran into while using KIDA for our research project was the inability of the code to accept a thermal rate coefficient function with more than three parameters (what we refer to as a "nonstandard" fitting function). Additionally, there was no clear way for the code to be run in a loop with auto-updating parameters, so we could understand what was happening over a large parameter space.

That's why we wrote pykida.

Look to the quickstart.py file to see what pykida can do. In it, we've defined a nonstandard thermal rate coefficient (in fact, this is the rate coefficient for the reaction of O on H3+ measured by de Ruette et al., which has been submitted to ApJ). We are able to run KIDA and Nahoon over a grid of temperatures from 10 to 400 K, and pipe the rate coefficient from Python into Fortran. We are able to do this with a small little trick: performing the initial computation in Python, and taking the value and storing it as the first coefficient of the Arrhenius-Kooij equation for that reaction, which then has the second and third coefficients set to 0.

This code will hopefully allow researchers to take advantage of more accurate fitting functions for their laboratory data, and it will allow for more efficient and automated astrochemical analysis of large parameter spaces. 