<img src="https://github.com/pekkosk/hotbit/blob/master/hotbit/doc/hotbit_logo_small.png" alt="logo" width="100">

# Hotbit 
Hotbit is an ASE density-functional tight-binding calculator that aims to provide
* an open-source DFTB code
* a handy companion for DFT (for easy & fast electronic structure analysis, for quick access to dynamical properties for testing, and for playing around)
* a compact and accessible code for everyone to inspect and modify (avoiding parallelization implies that the code is less suitable for large systems)
* an intuitive user interface (ideal for learning and teaching realistic electronic structure simulations)
* DFTB parametrization suite including interface to _libxc_ (see further instructions below)

  
## Take a closer look:

 1. [About hotbit](https://github.com/pekkosk/hotbit/wiki/About-hotbit)
 2. [Download and installation](https://github.com/pekkosk/hotbit/wiki/Download-and-installation)
 3. Manual:
   * [Calculator](https://github.com/pekkosk/hotbit/wiki/Calculator) (what Hotbit features)
   * [Getting started](https://github.com/pekkosk/hotbit/wiki/Getting-started)
   * [Charge self-consistency and Coulomb interaction](https://github.com/pekkosk/hotbit/wiki/Charge-self-consistency-and-Coulomb-interaction)
   * [Parameters and Parametrization](https://github.com/pekkosk/hotbit/wiki/Parameters-and-parametrization) (how to make parametrizations)
 4. [Code development](https://github.com/pekkosk/hotbit/wiki/Code-development) (for code developers)


If you find hotbit useful in your work, please cite ([pdf](http://users.jyu.fi/~pekkosk/resources/pdf/koskinen_CMS_09.pdf)):
```
    @article{koskinen_CMS_09,
      Author = {P. Koskinen, V. Mäkinen},
      Journal = {Computational Material Science},
      Title = {Density-functional tight-binding for beginners},
      Volume = {47},
      Pages = {237},
      Year = {2009}
    }
```

## environment variables
When installing with
```
python setup.py install --home=.
```
you can set the necessary environment variables by calling
```
bash env_exports
```

## libxc interface:
_libxc_ can be found at [www.tddft.org/programs/libxc](https://www.tddft.org/programs/libxc)
* For using the Hotbit Slater-Koster parametrization suite together with exchange-correlation functionals from _libxc_, you should install _libxc_ and its python module. For details on the installation, please refer to the instructions given with the _libxc_ package.
Currently, Hotbit supports LDA and GGA functionals (meta-GGAs and hybrids pending)
Specification of functionals is supported via their generic names. For a list of available functionals, please refer to
```
hotbit.parametrization.pylibxc_functionals.py
```
Support of direct specification via _libxc_ identifiers will be added soon.
