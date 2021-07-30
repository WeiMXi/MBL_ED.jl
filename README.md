# MBL_ED

[![Build Status](https://github.com/WeiMXi/MBL_ED.jl/workflows/CI/badge.svg)](https://github.com/WeiMXi/MBL_ED.jl/actions)
[![Coverage](https://codecov.io/gh/WeiMXi/MBL_ED.jl/branch/CPU/graph/badge.svg)](https://codecov.io/gh/WeiMXi/MBL_ED.jl)

It's a Julia program for repetition of the FIG. 2. in [Phys. Rev. Lett. 119, 075702 (2017) - Two Universality Classes for the Many-Body Localization Transition](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.075702)

## For start, just

```Julia
]add https://github.com/WeiMXi/MBL_ED.jl
using MBL_ED
# maybe you need install some pkg
# ]add JLD, Gadfly, DataFrames, Cairo, Statistics
# copy the example folder of the repository in the current directory
# then we can use the run.jl and data_process.jl
include("example/run.jl")   # run in parallel computing
include("example/data_process.jl")  # process data and get pdf of the result
```

Then, we can obtain result :

![result](./example/fin_MBL_CPU.svg)

## For adjust more parameters, change `example/run.jl`
