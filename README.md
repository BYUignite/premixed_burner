Burner stabilized premixed flame.

Corresponds to ISF 2a

Modified from cantera sample flamespeed.cpp. 

Needs Cantera with soot implementation (unless you turn off soot in burner.cc).

### Build
```
cd build
cmake ..
make
make install
```
* creates burner.x in the run directory


### Run in top level:
```
run/burner.x
```

### Outut
Creates burner.out.

### Other
* zT.dat is the fixed temperature profile that is interpolated
    * first column is z (m), second column is T (K)
* Experimental data: soot.exp, which has three experiments. 
    * These are separated in soot.exp_1, etc.



