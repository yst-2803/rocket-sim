# 1D Model Rocket Flight Simulator

A 1D flight simulator for an Estes Alpha III model rocket on a C6-7 motor, built from first principles in Python. It models altitude, velocity and apogee, and its predicted apogee agrees with OpenRocket to within 3.1%. That is a comparison between two simulations using the same drag coefficient, not a flight test. This project was inspired by my work experience at MBDA, which got me interested in how flight vehicles are simulated before they ever fly.

![Simulation vs OpenRocket](validation_vs_openrocket.png)

---

## The physics

The rocket is modelled as a point mass moving vertically, with upward taken as positive:

```math
m\frac{dv}{dt} = T(t) - mg - \tfrac{1}{2}\rho C_d A\, v|v|
```

where $T(t)$ is the measured thrust, $\rho$ is air density, $C_d$ is the drag coefficient and $A$ is the body's frontal area.

**Variable mass.** Applying conservation of momentum to the rocket and its exhaust gives

```math
m\frac{dv}{dt} = \dot{m}_e v_e - mg - D
```

so the thrust term is $T = \dot{m}_e v_e$.

$F = \frac{d(mv)}{dt}$ only holds for a system whose mass is constant, so it can't be applied to the rocket on its own: doing so adds an extra $v\,\dot{m}$ term whose size depends on the frame $v$ is measured in. The derivation therefore treats the rocket and its exhaust together as one closed system, and the thrust term $T = v_e\frac{dm_e}{dt}$ comes out of that momentum balance.

**Mass model.** Assuming constant exhaust velocity, $v_e = I_{total}/m_{prop}$, which gives

```math
m(t) = m_0 - m_{prop}\,\frac{I(t)}{I_{total}}, \qquad I(t) = \int_0^t T\,dt'
```

This means the mass at any time is the initial mass minus the propellant mass multiplied by the ratio of the impulse delivered so far to the total impulse.

**Drag direction.** Drag is written as $v|v|$ rather than $v^2$.

$v^2$ fixes the sign no matter what phase the rocket is in, as a square is always positive. Multiplying the velocity by its absolute value keeps its sign, so drag always opposes the motion whether the rocket is going up or coming down.

**Check:** with $g = 0$ and no drag, the equation reduces to the Tsiolkovsky rocket equation,

```math
\Delta v = v_e \ln\left(\frac{m_0}{m_f}\right)
```

The full derivation is in my handwritten notes: [derivation_notes.pdf](derivation_notes.pdf).

---

## Inputs

| Parameter | Value | Source |
|---|---|---|
| Motor | Estes C6-7 | [thrustcurve.org](https://www.thrustcurve.org) (certified test data) |
| Total impulse | 8.82 Ns | Thrust curve, trapezium rule |
| Launch mass $m_0$ | 57.1 g | 34 g rocket (Estes spec) + 23.1 g motor |
| Propellant mass $m_{prop}$ | 10.8 g | OpenRocket motor database |
| Dry mass (after burnout) | 46.3 g | 34 g rocket + 12.3 g empty casing |
| Body diameter | 25 mm ($A = 4.91 \times 10^{-4}$ m²) | Estes spec (the OpenRocket design uses its 24.79 mm BT-50 tube) |
| Drag coefficient $C_d$ | 0.56 | OpenRocket (rounded 1.5 mm fins) |
| Air density $\rho$ | 1.225 kg/m³ | Sea level, held constant (an ISA run is included for comparison) |
| Time step $\Delta t$ (Euler) | 0.01 s | — |

Conditions: no wind, vertical launch, sea level.

---

## Numerical methods

The state vector and its derivatives are

```math
\mathbf{y} = \begin{pmatrix} h \\ v \\ m \end{pmatrix}, \qquad
\frac{d\mathbf{y}}{dt} = \begin{pmatrix} v \\ \dfrac{T - mg - \tfrac{1}{2}\rho C_d A\, v|v|}{m} \\ -\dfrac{m_{prop}}{I_{total}}\,T \end{pmatrix}
```

Both solvers call the same derivative function, so any difference between them is purely numerical.

- **Euler (first order):** $\mathbf{y}_{n+1} = \mathbf{y}_n + \Delta t\,\mathbf{f}(t_n, \mathbf{y}_n)$
- **RK45:** SciPy's `solve_ivp` (Dormand–Prince, adaptive step), with `max_step = 0.01` s so it can't step over the thrust spike or burnout.
- **Thrust** between data points is found by linear interpolation.
- **Launch pad:** acceleration is held at zero while thrust < weight, so the rocket can't sink into the ground.
- **Apogee** is detected with an event function that finds where $v$ crosses zero.

The apogee event (which watches for $v$ crossing zero) fired at $t = 0$, because the rocket sits at exactly $v = 0$ on the pad until thrust exceeds its weight, which the solver counted as a crossing. I fixed this by only arming the event once the rocket is above 1 m.

---

## Verification: am I solving the equations correctly?

| Check | Result |
|---|---|
| Total impulse (trapezium rule) vs motor spec | 8.817 vs 8.8 Ns (0.2%) |
| Specific impulse $I_{total}/(m_{prop}\,g)$ from mass model vs spec | 83.2 vs ~83 s |
| Mass at burnout vs expected dry mass | 46.3 g vs 46.3 g (true by construction, so this only checks the mass integration) |
| Python Euler vs Excel Euler (same method, built independently) | 373.104 m vs 373.104 m |
| Euler ($\Delta t = 0.01$ s) vs RK45 | 373.10 vs 373.45 m (−0.09%) |

It tells us the time step used in Euler was small enough that the numerical error is negligible (about 0.1%). This means the gap to OpenRocket comes from the model's inputs and physics assumptions rather than the solver.

---

## Comparison with OpenRocket: am I solving the right equations?

OpenRocket is another simulation, and my drag coefficient is taken from it, so this tests the equations of motion and the integration rather than the aerodynamics. Proper validation needs flight data.

| Source | Apogee | vs OpenRocket |
|---|---|---|
| OpenRocket 24.12 | 385.5 m | — |
| **This simulator (RK45)** | **373.4 m** | **−3.1%** |
| This simulator (Euler) | 373.1 m | −3.2% |
| This simulator (RK45, ISA air density) | 376.8 m | −2.3% |
| This simulator (RK45, ISA air density and OpenRocket's body diameter) | 380.2 m | −1.4% |
| Analytical coast solution* | 378 m | −2.1% |
| Estes (manufacturer) | 351 m | −9.0% |

\*Exact coast-phase solution using OpenRocket's burnout state ($v_b = 117$ m/s, $h_b = 130$ m):

```math
\Delta h = \frac{m}{2k}\ln\left(1 + \frac{k v_b^2}{mg}\right), \qquad k = \tfrac{1}{2}\rho C_d A
```

| | This sim | OpenRocket | Difference |
|---|---|---|---|
| Max velocity | 116.2 m/s | 117 m/s | −0.7% |
| Time to apogee | 7.94 s | 8.12 s | −2.2% |

My simulation underestimates OpenRocket's apogee by 3.1%. I first put this down to constant air density and a constant drag coefficient in the coast. Breaking the gap down shows that is only part of it:

- **Air density explains about a quarter.** Switching to ISA density (1.225 kg/m³ on the pad, about 1.18 kg/m³ at apogee) raises the apogee from 373.4 m to 376.8 m, which closes 0.9 of the 3.1 percentage points.
- **Body diameter explains about another quarter.** I used the nominal 25 mm. The BT-50 body tube in the OpenRocket design (`AlphaIII_C6.ork`) is 24.79 mm, a 1.7% smaller frontal area. Using that diameter as well gives 380.2 m, closing another 0.9 points.
- **The rest is the drag coefficient changing with speed.** Working backwards from OpenRocket's own trajectory (`openrocket_effective_cd.py`) gives an effective $C_d$ between about 0.53 and 0.56 through the flight, against the constant 0.56 I used. With $C_d = 0.545$, ISA density and OpenRocket's diameter the apogee is 385.7 m (+0.1%). That value was picked after seeing the answer, so it shows where the gap comes from; it is not an independent prediction.
- **The gap opens during the burn, not only in the coast.** At burnout my simulation is already 1.1 m/s slower and 0.8 m lower than OpenRocket, at heights where air density has barely changed. The altitude plot above hides this; the residual plot below shows it.

![Residuals vs OpenRocket](residuals_vs_openrocket.png)

---

## Key findings

- In OpenRocket, changing the fins from square to rounded edges increased its predicted apogee by 43% ($C_d$ 1.08 → 0.56), so the drag coefficient is the most uncertain input to this model.
- My predicted apogee is sensitive to that input: changing $C_d$ by 0.03 moves it by about 11 m (3%), as much as the whole gap to OpenRocket.
- The peak thrust-to-weight ratio was about 25:1.
- At top speed, drag was about 5 times the rocket's weight (2.27 N against 0.45 N), so the total deceleration just after burnout is about 6 g.

---

## Limitations

- **Constant air density** in the baseline run (OpenRocket uses the International Standard Atmosphere; the ISA run above shows the effect).
- **Constant $C_d$**, taken from OpenRocket, so the comparison isn't fully independent: it tests the equations of motion and integration, not the aerodynamics.
- **1D only:** no wind, no weathercocking, no launch rod.
- **Constant exhaust velocity** assumed in the mass model.
- **No parachute:** after apogee the sim falls ballistically, so only the ascent is compared.

---

## Next steps

- Use OpenRocket's exported $C_d$ against speed instead of a constant
- Uncertainty analysis on $C_d$, mass and motor impulse
- Compare with real flight data from a logging altimeter
- Parachute descent and landing speed
- 2D trajectory with wind

---

## Running it

```bash
pip install -r requirements.txt
python rocket_sim.py
python openrocket_effective_cd.py
```

| File | Contents |
|---|---|
| `rocket_sim.py` | Simulation, solvers and plots |
| `Estes_C6.csv` | Certified C6 thrust curve (thrustcurve.org) |
| `openrocket_effective_cd.py` | Works out the drag coefficient OpenRocket was effectively using |
| `openrocket_C6-7.csv` | OpenRocket export used for the comparison |
| `AlphaIII_C6.ork` | OpenRocket design file the export came from |
| `validation_vs_openrocket.png` | Comparison plot |
| `residuals_vs_openrocket.png` | Difference between this simulator and OpenRocket against time |
| `derivation_notes.pdf` | Handwritten derivations |
| `excel_verification.xlsx` | Hand verification: impulse, mass model, coast estimate and Excel Euler |
