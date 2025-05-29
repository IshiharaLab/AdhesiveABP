# Attractive ABP

Attractive active Brownian particle model used in S. Shimamura et al. [1]


---

## Description

The source code for simulating the model. This repository contains two main simulation codes and a library:
- `LJ_ABP.cpp`: Standard simulation for attractive ABPs in a square box.
- `LJ_ABP_rectangle.cpp`: Simulation in a rectangular box, used for studying cluster escape dynamics.
- `lib`: a library for constructing the cell list and visualizing results with Gnuplot.

---


## Requirements

* g++ with C++17 support  
* GNU Scientific Library (GSL) development headers (`libgsl-dev`)
* Gnuplot (for visualization, optional)

---

## Usage

1. Clone the repository and inspect the directory structure:  
   ```text
   project/
   ├── LJ_ABP.cpp
   ├── LJ_ABP_rectangle.cpp
   ├── lib/
   │   ├── particle_neighbours.cpp
   │   └── adhABPlib.c
   └── README.md
2. Set the parameters at the top of `LJ_ABP.cpp` (or `LJ_ABP_rectangle.cpp`). See `Parameter setting` to reproduce the specific data for the results in [1].
3. Use the following command to compile the program (if you use `LJ_ABP_rectangle.cpp`, replace `LJ_ABP.cpp` with it):
```sh
g++ -std=c++17 LJ_ABP.cpp lib/particle_neighbours.cpp lib/adhABPlib.c \
    -lgsl -lgslcblas -lm -o LJ_ABP
```

4. Run the simulation and redirect output to a file:
```sh
mkdir -p results
./LJ_ABP <f0> <u0> <seed> > results/output.txt
<f0>: active force magnitude
<u0>: adhesion strength
<seed>: random number seed
```

5. To produce an animated GIF of the particle trajectories: Edit GNUPLOT_SAVE_GIF to 1 in `LJ_ABP.cpp` (or `LJ_ABP_rectangle.cpp`).
   Recompile as above and run the simulation; then `testanime.gif` will be generated in the working directory.



---
## Parameter setting
Unless explicitly stated, all other parameters are fixed as in the source code. Below is the mapping of figures/movies to code files and parameter settings.

---

### 1. Most part of the results  
**Figures**: Fig.1-(b),(c), Fig.2-(a),(b), Fig.3-(b), Fig.S1-(a),(c), Fig.S2-(a),(b),(d), Fig.S3-(a),(b),(d), Fig.S5-(a),(b)  
**Movies**: Movie-1, 2, 3, 4  
**Code**: `LJ_ABP.cpp`  
**Parameters**:
```
INITIAL_CONDITION_TYPE = 1
(N, LX, LY) = (5000, 100, 100)
SIMULATION_TIME = 40000
seed = 1
f0 = 0.15, 0.3, ..., 2.55, 2.7
u0 = 0, 2.5, ..., 15.0
```

---

### 2. For measuring cluster activity  
**Figures**: Fig.3-(a), Fig.S2-(c), Fig.S3-(c)  
**Code**: `LJ_ABP.cpp`  
**Parameters**:
```
INITIAL_CONDITION_TYPE = 1
SIMULATION_TIME = 100000
seed = 2, 3, ...,9
(f0, u0) = (0.3, 0.0), (1.2, 0.0), (1.2, 10.0)
```

---

### 3. For visualizing cluster escape  
**Figures**: Fig.3-(c),(d)  
**Movies**: Movie-5, 6  
**Code**: `LJ_ABP_rectangle.cpp`  
**Parameters**:
```
INITIAL_CONDITION_TYPE = 1
(N, LX, LY) = (6250, 250, 50)
SIMULATION_TIME = 40000
seed = 1
(f0, u0) = (1.2, 0.0), (1.2, 10.0)
```

---

### 4. For large system size study  
**Figures**: Fig.S1-(b), Fig.S4-(b)  
**Code**: `LJ_ABP.cpp`  
**Parameters**:
```
INITIAL_CONDITION_TYPE = 1
SIMULATION_TIME = 40000
seed = 1
(N, LX, LY) = (5000, 100, 100), (11250, 150, 150), (20000, 200, 200)
Fig.S1-(b): f0 = 0.1, 0.2, ..., 3.0, u0 = 10
Fig.S4-(b): (f0, u0) = (0.3, 0.0), (1.2, 0.0), (1.2, 10.0)
```

---

### 5. For assessing AG-MC crossover  
**Figure**: Fig.S4-(a)  
**Code**: `LJ_ABP.cpp`  
**Parameters**:
```
INITIAL_CONDITION_TYPE = 0
SIMULATION_TIME = 40000
seed = 1
(N, LX, LY) = (5000, 100, 100)
f0 = 0.6
u0 = 0.0, 1.0, ..., 6.0
```

---

## References
 [1] Sota Shimamura, Nen Saito, Shuji Ishihara. Attraction-Induced Cluster Fragmentation and Local Alignment in Active Particle Systems. doi.org/10.48550/arXiv.2505.19118



```python

```
