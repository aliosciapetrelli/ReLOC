# ReLOC
## Pairwise registration by local orientation cues

[ ![License] [license-image] ] [license]

[license-image]: https://img.shields.io/badge/license-gpl-green.svg?style=flat
[license]: https://github.com/aliosciapetrelli/ReLOC/blob/master/LICENSE

Description
-----------
C++ implementation of the pairwise 3D registration algorithm proposed in:

[Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015.](http://onlinelibrary.wiley.com/doi/10.1111/cgf.12732/epdf)

Webpage
-----------
http://www.vision.deis.unibo.it/research/78-cvlab/109-reloc

Usage
-----------
An example can be found in [ReLOC_testmain.cpp](https://github.com/aliosciapetrelli/ReLOC/blob/master/ReLOC_testmain.cpp).

The algorithm can be evaluated through the benchmark proposed in the same paper. The C++ code of the benchmark can be downloaded from [here](https://github.com/aliosciapetrelli/Pairwise3DRegistrationEvaluation) and applied to ReLOC as explained in [ReLOC_benchmarkmain.cpp](https://github.com/aliosciapetrelli/ReLOC/blob/master/ReLOC_benchmarkmain.cpp).

Dependencies
-----------
The framework requires the [VTK](http://www.vtk.org/) library.

The code has been tested with VTK 5.10 on Windows 7 and Microsoft Visual Studio 2010.
