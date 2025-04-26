## Optimizing damper placement in automotive mechanical suspension using Bayesian Optimization

### Overview
* Parameters of an automotive suspension system are optimized to minimize an objective function that balaances:
  * Comfort (Reducing acceleration felt by passengers)
  * Vibration (Minimizing chassis displacement)
  * Handling (Reducing pitch, roll, and tire force variations)
  * Constraints (Ensuring suspension stroke stays within limits)
* Damper placement is indirectly optimized through parameters like motion ratio, inclination angle, and damping coefficients which affect how damping forces are transmitted to the suspension system.

### Functions
* Simulates the half-car model and evaluates the objective function for a given set of suspension parameters.
* Compression and rebound damping coefficients (Ns/m) controls damper force.
* Motion ratio scales the damper's effect based on its mechanical advantage in the suspension linkage.
* Inclination angle affects the effective damping force via cosine projection while knee point, blowoff, hysteresis, temp coeff. define nonlinear damping behavior.
* Modelling the vehicle's response to road inputs using a 10-state half-car model (sprung/unsprung masses, pitch, roll, etc.).
* Damping forces are computed incorporating nonlinear effects like blow-off thresholds and temperature-dependent damping.
* Outputs (acceleration, displacement, pitch, roll, tire forces, max stroke) are used to compute the objective function
* Damper placement affects comfort (via acceleration) and handling (via pitch/roll).
* Bayesian Optimization is utilized and it optimizes the suspension parameters by modeling the objective function with a Gaussian Process (GP) and maximizing the Expected Improvement (EI) acquisition function.
* Initialization uses Latin Hypercube Sampling to generate initial parameter sets and optimization loop fits the GP to observed data, predicts mean/variance for new points, and selects the next parameters by maximizing EI over 1000 random samples.
* All 16 parameters are optimized within bounds.
* The GP indirectly learns how motion ratio and inclination angle impact the objective, optimizing their values to improve damping effectiveness.
* A random road profile is generated to simulate realistic road inputs.
* An ISO 8608-like roughness model with smoothed noise, providing displacement and velocity inputs to the half-car model is used.
* The road profile drives the suspension dynamics, testing the damper’s ability to absorb disturbances.
* Power Spectral Density (PSD) of chassis displacement to evaluate vibration, which is influenced by damper settings is computed.

### Optimization of damper placement
* Motion ratio (0.5–1.5) determines how damper motion translates to suspension motion. A higher ratio increases damping force but may reduce stroke range.
* Inclination angle (0–30°) affects the effective damping force. A smaller angle (closer to vertical) maximizes force transmission, while a larger angle may reduce wear or fit packaging constraints.
* Damping coefficients control the damper’s force-velocity relationship, tuned to balance comfort and handling.
* Non-linear damping parameters i.e blowoff, knee point, hysteresis, temp coefficient allow the damper to adapt to different velocities and conditions, indirectly influenced by placement.
* Bayesian optimizer explores combinations of these parameters, using the GP to model the objective function’s response.
* The EI acquisition function prioritizes parameter sets that are likely to improve the objective, efficiently navigating the trade-offs between comfort, vibration, handling, and constraints.

### Sample workflow
* 10 initial parameter sets are sampled (e.g., $C_c$ = 2000, $C_r$ = 4000, motion ratio = 1.0, inclination = 15°).
* The `VehicleDynamics::evaluateObjective` function simulates each set, computing the objective value (e.g., 5.2).
* The GP is fitted to the initial (parameters, objective) pairs.
* For each iteration (30 total):
  * 1000 random parameter sets are evaluated using the EI function.
  * The set with the highest EI (e.g., $C_c$ = 2500, $C_r$ = 4500, motion ratio = 1.2, inclination = 10°) is simulated.
  * The new result updates the GP, refining predictions.
* After 30 iterations, the best parameters are returned (e.g., $C_c$ = 3000, $C_r$ = 5000, motion ratio = 1.1, inclination = 8°, objective = 4.8).
* These values optimize damper placement (via motion ratio and inclination) and damping behavior for the given road profile and vehicle model. 
