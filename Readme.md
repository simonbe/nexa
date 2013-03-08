Nexa is an experimental neural simulator for parallel simulation of large-scale neural network models at a high level of biological abstraction and for exploration of the simulation methods involved. 

It includes among other things:
- Parallelization using the Message Passing Interface (MPI).
- Firing-rate models.
- Capabilities to build networks using machine learning inspired methods for e.g.  self-organization of network architecture and for structural plasticity.
- Analysis running in parallel beside the simulation.


Also see

Benjaminsson, S. and Lansner, A. (2012). [Nexa: A scalable neural simulator with integrated analysis] (http://informahealthcare.com/doi/abs/10.3109/0954898X.2012.737087). Network: Computation in Neural Systems, 23(4): 254-271.


### Installation ###

[Download the latest version as a zip file](https://github.com/simonbe/nexa/zipball/master) or clone the repo:

```bash
$ git clone git@github.com:simonbe/nexa
```

For easy access on Windows, [GitHub for Windows](http://windows.github.com/) can be used.

**Project files and Makefiles**

Repo includes a Visual Studio 2010/2012 project file. To get MPI libraries and parallel debugging, [Microsoft Compute Cluster Pack SDK](http://www.microsoft.com/en-us/download/details.aspx?id=239) or [Microsoft HPC Pack](http://www.microsoft.com/en-us/download/details.aspx?id=8433) are recomended.
For other platforms (currently CRAY XE6 using the GNU compiler and BG/L and BG/P using the IBM compiler), Makefiles exist in the [Resources/Makefiles](https://github.com/simonbe/nexa/tree/master/Resources/Makefiles) directory.

### Setting up a network ###

[Step-by-step network setup and simulation](https://github.com/simonbe/nexa/wiki/Step-by-step-network-setup-and-simulation)

### Modifying and adding models, e.g. for synaptic plasticity ###

[Adding a model](https://github.com/simonbe/nexa/wiki/Adding-a-model)

### License ###
Nexa is licensed under the [GNU LGPL license](http://www.gnu.org/licenses/lgpl.html).
