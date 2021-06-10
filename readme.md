# Thermodynamic library: humidityRhoThermo
This library takes humidity effects into account and is mainly built for HVAC analysis. The validity is limited to around 1 bar and -50 degC to 100 degC. The library is built for different OpenFOAM versions. For more complex analysis (different pressure than atmospheric one or higher temperatures) one should use the reactingFoam for taking into account the humidity.

# How to use it

Clone the repository to any place you want using the following command:
```console
@-: git clone https://github.com/shor-ty/humidityRhoThermo.git
```

After that load your OpenFOAM environment (if not already happend) and move into the repository. Here checkout your version you want:
```console
@~: git checkout <TAB><TAB>
OpenFOAM-6.x
OpenFOAM-7.x
OpenFOAM-v1712
OpenFOAM-v1806
OpenFOAM-v1812
OpenFOAM-v1906
OpenFOAM-v1912
OpenFOAM-v8
@~: git checkout OpenFOAM-v8
```
After you switched to your OpenFOAM version, you have to compile the library first and after that the solver:
```console
@~: cd src/thermodynamic/basic/
@~: wmake libso
@~: cd ../../../applications/solvers/heatTransfer/buoyantHumidityPimpleFoam/
@~: wmake
```
You are done. After that you can use the buoyantHumidityPimpleFoam solver.

# Sponsored
This project was sponsered by Tian Building Engineering
