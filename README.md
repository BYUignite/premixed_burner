Burner stabilized premixed flame.

Corresponds to ISF 2a

Modified from cantera sample flamespeed.cpp. 

### Build
```
cd build
cmake ..
make
make install
```
* creates burner.x in the top directory


### Run in top level:
```
./burner.x
```

### Output:
* burner.out

### Other
* zT.dat is the fixed temperature profile that is interpolated
    * first column is z (m), second column is T (K)



