# Parameter setting

This is the parameter setting to reproduce the specific data for the results in [1] and its Supplemental Material. See `README.md` for the usage for the source codes.  

Unless explicitly stated, all other parameters are fixed as in the source code. Below is the mapping of figures/movies to code files and parameter settings.

---

## 1. Most part of the results  
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

## 2. For measuring cluster activity  
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

## 3. For visualizing cluster escape  
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

## 4. For large system size study  
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

## 5. For assessing AG-MC crossover  
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
