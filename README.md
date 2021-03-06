# Gray-Scott Reaction-Diffusion

An implementation in Processing by Kaylah Facey.

## Commands:

* i - Initialize the system with a fixed rectangular region that has specific u and v concentrations.
* space bar - Start or stop the simulation (toggle between these).
* u - At each timestep, draw values for u at each cell (default).
* v - At each timestep, draw values for v at each cell (rather than u).
* d - Toggle between performing diffusion alone or reaction-diffusion (reaction-diffusion is default).
* p/0 - Toggle between constant f, k for each cell (default) and spatially-varying f, k.
* 1 - Set parameters for spots (f = 0.035, k = 0.0625)
* 2 - Set parameters for stripes (f = 0.035, k = 0.06)
* 3 - Set parameters for spiral waves (f = 0.0118, k = 0.0475)
* 4-7 - Set parameters for custom patterns:
-   4 - concentric squares (f = .0271, k = .0557)
-   5 - flower (f = .0624, k = .0614)
-   6 - blinking spots (f = .0140, k = .0503)
-   7 - disappearing spots (f = .0108, k = .0531)
* (mouse click) - Print values for u and v at cell.  If in spatially-varying parameter mode, also print values for k and f at the cell.

# Custom Commands:
* t - Set boundary condition to periodic (toroidal, default).
* z - Set boundary condition to dirichlet (fixed-value).
* f - Set boundary condition to neumann (zero-derivative).
* >/. - Increase dt by .1 (up to 3; default 1)
* </, - Decrease dt by .1 (down to .1; default 1)
