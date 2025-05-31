# Initial Conditions

The ```GeoModBox.jl``` contains several [routines](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/InitialCondition/2Dini.jl) or [structures](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Structures.jl) to define certain parameters or to setup a certain initial anomaly, either for properties defined on their correspondig grid (i.e., temperature, velocity, or phase) or for tracers. Within the [examples] and the [exercise] one can choose different initial temperature and velocity conditions.  

## Mutable Structures

Implementing the scaled version of the thermal convection into the ```GeoModBox.jl``` showed that immutable ```NamedTuples```are numerically unfavourable for the scripts. Some initially defined parameters, like the model height, need to be modified in order to scale them. 

```Geometry()```

```Physics()```

```GridSpacing()```

```DataFields()```

```TimeParameter```

### Initial Temperature

```IniTemperature!()``` 

### Initial Velocity

```IniVelocity!()``` 

### Initial Phase

```IniPhase!()```

## Initialize Tracers

