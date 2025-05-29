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
## Parameters

 See Paremeter_setting.md for parameter values reproducing the results described in  S. Shimamura et al. [1].


---

## References
 [1] Sota Shimamura, Nen Saito, Shuji Ishihara. Attraction-Induced Cluster Fragmentation and Local Alignment in Active Particle Systems. doi.org/10.48550/arXiv.2505.19118



```python

```
