# Light Forward Ray Tracing in Kerr Spacetime

This code calculates the trajectory of light (ray - geometrical optics) on a background Kerr spacetime


The code solves a set of ODEs numerically. These equations are based on the original works of [Mathisson 1937](http://inspirehep.net/record/48323/citations), [Papetrou 1951](https://royalsocietypublishing.org/doi/10.1098/rspa.1951.0200) and [Dixon 1964](https://link.springer.com/article/10.1007%2FBF02734579). Consequently these equations are known as the MPD equations. 
More recent works can be found in [Mashoon & Singh 2006](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.74.124006), [Singh 2005](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.72.084033) and [Singh, Wu & Sarty 2014](https://academic.oup.com/mnras/article/441/1/800/983192).

Additional interesting discussion on the motion of extended bodies in GR can be found in [Consta and Natatio, 2015](https://arxiv.org/abs/1410.6443)


## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This code is written in FORTRAN with a [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. **Other compilers have not been tested.** The gfortran installation binaries can be found [here](https://gcc.gnu.org/wiki/GFortranBinariels), although typically gfortran comes pre-installed on most Linux/Unix systems. If you have [Homebew](https://brew.sh/) installed on OSX, you can simply run 


```
brew install gcc
```

### Starting steps
After [cloning the repo](https://help.github.com/en/articles/cloning-a-repository), the first thing to do is to set the path to the output files that the code will produce.
This can be done by first setting the path destination in `config.sh' and then running the executable. Note taht sometimes permissions on `.sh' files will need to be changes as,

echo 'export RayTracingDir="/Users/tomkimpson/Data/RayTracing3/"' >> ~/.bash_profile
source ~/.bash _ profile



Set this to point to a local direcory.


check it is there with env or vim ~/.bashprofile


The code should then run as is, out of the box. Try

```
run.py
```

to compile and run the code. Once you have checked that everything is running OK, you can then start playing. The code structure (mdoules, subroutines etc.) is outlined below.


If making edits to the code, try to keep to the [FORTRAN Style Guide](https://www.fortran90.org/src/best-practices.html)

## Structure

`parameters.f` defines all the system parameters. That is, anything that needs changing (e.g. eccentricity, orbital period, BH mass) can be modified in this module


`constants.f` is for calculations with those parameters for use later in the code. It can effectively be ignored - no changes should be necessary to this file

`main.f` is where the code program is run from. After setting up the initial conditions (`initial_conditions.f`) it then goes on to integrate the equations and save the output (`rk.f`). 

In turn, `rk.f` calls `derivatives.f` and `tensors.f` to calculate e.g. the curvature tensors, ODEs and then integrates numerically.


A python wrapper has been provided to compile and run the code, `run.py`. We use a `-O3` optimization. [See the docs](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) for discussion on the optimization flags


## Numerical Method
We integrate the equations using a Runge-Kutta-Fehlberg algorithm with adaptive stepsize. See [Press et al.](https://dl.acm.org/citation.cfm?id=141273)


### Accuracy tests
When integrating numerically, an important consideration is the accuracy of the method. We can assess this by independelty evaluating the conserved quantities $E,L,Q$ independently.



## ToDO()

* Add E/L/Q optinal checks
* Add GW module
* fix z outout of RK


## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc


