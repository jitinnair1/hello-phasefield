This is a collection of some simple phasefield problems solved in MATLAB

# Contents

0. Some basic exercises to get started with MATLAB
```
├── basic_exercises
│   ├── IdealSolution.m
│   ├── PhaseDiagram.m
│   ├── RegularSolution.m
│   ├── RegularSolutionFunction.m
│   ├── StirlingPlot.m
│   └── fft_basics.m
```

1. Solving 1D difusion equation using explicitlty using Finite Difference Methods (using different boundary conditions) and Fourier Transforms
```
├── 1D_diff_explicit
│   ├── DiffExplicitFDMv1.m
│   ├── DiffExplicitFDMv2.m
│   ├── DiffExplicitFDMv3.m
│   └── DiffExplicitSpectral.m
```

2. Solving 1D difusion equation using implicitly using Finite Difference Methods (using different boundary conditions) and Fourier Transforms
```
├── 1D_diff_implicit
│   ├── DiffImplicitFDMv1.m
│   ├── DiffImplicitFDMv2.m
│   ├── DiffImplicitFDMv3.m
│   └── DiffImplicitSpectral.m
```

3. Solving the Allen Cahn equation using Fourier Transforms
```
├── allen_cahn
│   ├── AC_1D_Spectral.m
│   ├── AC_2D_Spectral.m
│   └── GibbsThomson.m
```

4. Solving the Cahn-Hilliard equation using Fourier Transforms
```
├── cahn_hilliard
│   ├── CH_1D_FDM.m
│   ├── CH_1D_Spectral.m
│   ├── CH_2D.m
│   ├── CH_FDM_GridPoint.m
│   ├── CH_Photo.m
│   ├── CH_SpectalvsFDM.m
│   ├── CH_Variation_Analysis_A.m
│   ├── CH_Variation_Analysis_kappa.m
│   └── CH_intuition.m
```

5. Spinodal decomposition
```
└── spinodal
|   ├── Spinodal2D.m
|   ├── Spinodal2D_EnergyPlots.m
|   ├── Spinodal3D.m
|   ├── SpinodalPhaseDiagram.m
|   └── write_vtk_grid.m
```

6. Fan-Chen model for two grain system
```
├── fan_chen
│   └── FanChen.m
```

7. How the solution to the Diffusion equation is affected by choice of 'dx' and 'dt'
```
├── gridsize_effect
│   ├── DiffExplicinFunc.m
│   └── GridPointVariation.m
```

8. Calculation of interfacial width using Cahn-Hiliard equation
```
├── interfacial_width
│   ├── CH_Interfacial_FDM.m
│   ├── CH_Interfacial_Spectral.m
│   ├── CH_Interfacial_width_direct.m
│   ├── CH_Interfacial_width_func.m
│   ├── CH_Interfacial_width_slope_v0.m
│   ├── CH_Interfacial_width_slope_v1.m
│   ├── CH_Interfacial_width_slope_v2.m
│   ├── CH_Interfacial_width_slope_v3.m
│   └── CH_Interfacial_width_slope_variation.m
```

9. Solving sintering of two particle system
```
├── sintering
│   ├── free_energy_diff.m
│   ├── micro_sint_pre.m
│   ├── sintering.m
│   └── write_vtk_grid_values.m
```


# References
- [NPTEL course on Phasefield Modelling](https://nptel.ac.in/courses/113/101/113101072/)
- [Programming Phase-Field Modeling by S. Bulent Biner](https://www.springer.com/gp/book/9783319411941)
