=====================================
This is MeshKit, an open-source library for mesh generation functionality.

MeshKit is distributed under the terms of the GNU LGPL. See the LICENSE.txt

The example/ directory contains the examples from the described in
the doxygen page (http://www.mcs.anl.gov/~fathom/meshkit-docs/html/examples.html).
The test/ directory contains tests for individual algorithms and tools.

MeshKit user and developer guide are available online:
http://www.mcs.anl.gov/~fathom/meshkit-docs/html/userguide.html
http://www.mcs.anl.gov/~fathom/meshkit-docs/html/devguide.html

=====================================
Installation Instructions:
=====================================

1. System requirements
=====================================
MeshKit is currently supported on Linux and Linux-like operating systems (including MacOS); support for Microsoft Windows is under development and should be available by the next release.

2. Configure/Build/Installtion of MeshKit
=====================================
2.1 Prerequisites
=====================================
MeshKit requires the following libraries to be installed before configuration:

• CGM: a library for representation, query and modification of geometric
   models; see https://trac.mcs.anl.gov/projects/ITAPS/wiki/CgmFromScratch
   for details on obtaining and building CGM.
• MOAB: a library for representing structured and unstructured mesh; see 
   https://trac.mcs.anl.gov/projects/ITAPS/wiki/BuildingMoab
   for details on obtaining and building MOAB.
• Autotools: this is a set of Linux utilities for configuring software packages.
   Autotools can be found in most Linux package managers, and usually consists
   of the Autoconf and Automake packages.

In addition, if a parallel version of MeshKit is desired, the user must have the Message
Passing Interface (MPI) available on their computer; binary versions of MPI can be found in
most Linux package managers.

2.2 Download, Configure, Build, Install
=====================================
MeshKit source code is maintained in a world-readable svn repository, located at
https://svn.mcs.anl.gov/repos/fathom/MeshKit/trunk/. By default, MeshKit uses a GNU
Autotools-based configuration process. The following steps should be used to configure,
build, and install MeshKit:

• Unpack the source tarball into a directory referred to below as <MK_DIR> and
   change directory into that location.
• Execute ‘autoreconf –fi’. This executes a series of tools in the autotools suite, storing some generated files in the ‘config’ subdirectory.
• Execute ‘./configure’ with appropriate options. Three configure options are required, specifying the locations of CGM (--with-igeom=<location>), MOAB (--with-imesh=<location>) and LASSO (--with-irel=<location>). 
Other useful configure options are the installation location (--prefix=<location>) and specifying debug or optimized builds (--enable-debug, --enable-optimized, respectively). 
For a complete list of options, execute the command ‘./configure –help’. 
After a successful configuration, a set of Makefile’s are generated in the proper subdirectories.
• To complete the build of MeshKit, execute ‘make’.
• To install MeshKit, execute ‘make install’. 
If the install location was not specified on the configure line, one can specify a location in this step by using the command ‘make prefix=<location> install’.

For those wishing to use the Python interface, MeshKit and its dependencies should be
configured to build shared libraries, using the ‘--enable-shared’ configure option where
appropriate.
Once the MeshKit library has been built, it is ready for inclusion into user-developed
applications (any MeshKit-packaged programs, e.g. those that constitute RGG, will be
installed in the ‘bin’ directory). To aid in building user-developed applications, MeshKit also
writes a file ‘meshkit.make’, which can be included directly into application makefiles. This
file defines the following make variables useful for building MeshKit-based applications:
• MESHKIT_INCLUDES, MESHKIT_CPPFLAGS: compiler options
   pointing to all directories containing include files available to applications,
  including those for CGM and MOAB; also, CPP definitions controlling which
 optional external meshing tools have been configured into MeshKit.
• MESHKIT_LIBS_LINK: linker options necessary to satisfy all functions
   included in MeshKit.
The ‘examples’ subdirectory in the MeshKit source installation contains an example
makefile showing how these make variables can be used to compile and link MeshKit-based
applications.

3. Release MeshKit: Version 1.1
=====================================
This version has been tested on Ubuntu 10.04, 12.04 and Mac OS X 10.6.8 to 10.9.
Some dependencies may have restrictive license. 

3.1 Components
=====================================
The following version of dependencies are found to be compatible with MeshKit v 1.1
MPICH2 1.4.1p1, 1.5, open-mpi 1.6.4
MOAB Git SHAR1: bd52ba12517416f4b6d2162696a41583b73d52ed
CGM  https://svn.mcs.anl.gov/repos/ITAPS/cgm/tags/13.1.1
LASSO https://svn.mcs.anl.gov/repos/ITAPS/Lasso/trunk Rev: 6103

TRIANGLE: Triangle provides a 2D triangle mesh generation algorithm.

CAMAL/CUBIT 5.1.0: Sandia National Lab’s proprietary CAMAL library provides access to tetrahedral, tri-advance and paving algorithms.

NetGen 4.9.13: An open source automatic tetrahedral mesh generation library.

Mesquite: Open source mesh optimization package MESQUITE 

IPOPT: 3.10.0 (Coin), MA27 solver http://www.hsl.rl.ac.uk/download/coinhsl-archive/2011.10.03/  
Ipopt an optimization library used by interval assignment algorithm developed in MeshKit.


