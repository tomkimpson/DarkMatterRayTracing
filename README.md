# Light Ray Tracing in the a Dark Matter Spacetime

This code calculates the trajectory of light (ray - geometrical optics) on a background Kerr spacetime surrounded by a dark matter halo.


## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This code is written in FORTRAN with a [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. **Other compilers have not been tested.** The gfortran installation binaries can be found [here](https://gcc.gnu.org/wiki/GFortranBinariels), although typically gfortran comes pre-installed on most Linux/Unix systems. If you have [Homebew](https://brew.sh/) installed on OSX, you can simply run 


```
brew install gcc
```

### Starting steps
After [cloning the repo](https://help.github.com/en/articles/cloning-a-repository), the first thing to do is to set the path to the output files that the code will produce.
This can be done by setting the environment variable as

```
echo 'export DarkMatterDir="/Users/tomkimpson/Data/DM/"' >> ~/.bash_profile
source ~/.bash_profile
```

Just change the path `Users/tomkimpson/Data/DM/` to some appropriate local path. 

You can check the environemnt variable has been added to `bash_profile` by either `env` or `vim ~/.bashprofile`


The code should then run as is, out of the box. Try

```
run.py
```
to compile and run the code. Once you have checked that everything is running OK, you can then start playing. The code structure (modukes, subroutines etc.) is outlined below.


If making edits to the code, try to keep to the [FORTRAN Style Guide](https://www.fortran90.org/src/best-practices.html)

## Structure

### src
`parameters.f` defines all the system parameters. That is, anything that needs changing (e.g. BH mass, BH spin) can be modified in this module

`constants.f` is for calculations with those parameters for use later in the code. It can effectively be ignored - no changes should be necessary to this file

`main.f` is where the code program is run from. After setting up the initial conditions (`initial_conditions.f`) it then goes on to integrate the equations and save the output (`integration.f') 

`tensors.f` contains some useful subroutines for calculating e.g. metric, vector magnitudes etc.


A python wrapper has been provided to compile and run the code, `run.py`. We use a `-O3` optimization. [See the docs](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) for discussion on the optimization flags

### tools

`PlotRays.py` does what is says on the tin. Can be switched between 3d and 2d by changing the opening `d' parameter


## Numerical Method
We integrate the equations using a Runge-Kutta-Fehlberg algorithm with adaptive stepsize. See [Press et al.](https://dl.acm.org/citation.cfm?id=141273)



## Running in parallel 
The code is compiled using `-fopenmp` to allow for parallel processing. Geodesic ray tracing naturally lends itself very well to parallel computation. See the [OpenMP docs](https://www.openmp.org/wp-content/uploads/openmp-4.5.pdf). To set the number of threads, use

```
export OMP_NUM_THREADS=1
```


### Accuracy tests
When integrating numerically, an important consideration is the accuracy of the method. We can assess this by independently evaluating the Carter Constant, *Q*


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



