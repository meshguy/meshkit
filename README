# MeshKit

MeshKit is an open-source library for mesh generation functionality. MeshKit is distributed under the terms of the GNU LGPL. See the LICENSE.txt file in the repository.

The example/ directory contains the examples from the described in
the [Doxygen page]. The test/ directory contains tests for individual algorithms and tools.

MeshKit user and developer guide are available [online].

# Installation Instructions

## 1. System requirements

MeshKit is currently supported on Linux and Linux-like operating systems (including MacOS); support for Microsoft Windows is under development and should be available by the next release.

## 2. Configure/Build/Installtion of MeshKit

### 2.1 Prerequisites

MeshKit requires the following libraries to be installed before configuration:

  - CGM: a library for representation, query and modification of geometric models; see [CGM build instructions] for details on obtaining and building CGM.
  - MOAB: a library for representing structured and unstructured mesh; see [MOAB build instructions] for details on obtaining and building MOAB.
  - Lasso: a library for recovering and querying relations between mesh and geometry; see [Lasso build instructions] for details on obtaining and building Lasso.
  - Autotools: this is a set of Linux utilities for configuring software packages. Autotools can be found in most Linux package managers, and usually consists of the Autoconf and Automake packages.

In addition, if a parallel version of MeshKit is desired, the user must have the Message Passing Interface (MPI) available on their computer; binary versions of MPI can be found in most Linux package managers.

### 2.2 Download, Configure, Build, Install

MeshKit source code is maintained in a world-readable git repository, located at
https://bitbucket.org/fathomteam/meshkit.git. By default, MeshKit uses a GNU
Autotools-based configuration process. The following steps should be used to configure, build, and install MeshKit:

  - Unpack the source tarball into a directory referred to below as <MK_DIR> and change directory into that location.
  - Execute ‘autoreconf –fi’. This executes a series of tools in the autotools suite, storing some generated files in the ‘config’ subdirectory.
  - Execute ‘./configure’ with appropriate options. Three configure options are required, specifying the locations of CGM (--with-igeom=<location>), MOAB (--with-imesh=<location>) and LASSO (--with-irel=<location>).

Other useful configure options are the installation location (--prefix=<location>) and specifying debug or optimized builds (--enable-debug, --enable-optimized, respectively). 

For a complete list of options, execute the command ‘./configure –help’. 
After a successful configuration, a set of Makefile’s are generated in the proper subdirectories.

  - To complete the build of MeshKit, execute ‘make’.
  - To install MeshKit, execute ‘make install’. 

If the install location was not specified on the configure line, one can specify a location in this step by using the command ‘make prefix=<location> install’.

For those wishing to use the Python interface, MeshKit and its dependencies should be
configured to build shared libraries, using the ‘--enable-shared’ configure option where appropriate.

Once the MeshKit library has been built, it is ready for inclusion into user-developed applications (any MeshKit-packaged programs, e.g. those that constitute RGG, will be installed in the ‘bin’ directory). To aid in building user-developed applications, MeshKit also writes a file ‘meshkit.make’, which can be included directly into application makefiles. This file defines the following make variables useful for building MeshKit-based applications:
  
  - MESHKIT_INCLUDES, MESHKIT_CPPFLAGS: compiler options pointing to all directories containing include files available to applications, including those for CGM and MOAB; also, CPP definitions controlling which optional external meshing tools have been configured into MeshKit.
  - MESHKIT_LIBS_LINK: linker options necessary to satisfy all functions included in MeshKit.

The ‘examples’ subdirectory in the MeshKit source installation contains an example
makefile showing how these make variables can be used to compile and link MeshKit-based applications.

**Note: Automatic CMake-based build that builds MeshKit, MOAB, CGM, Lasso, OCE, FTGL, HDF5 can be found in folder <MeshKit>/trunk/doc/SuperBuild

[Doxygen page]: http://ftp.mcs.anl.gov/pub/fathom/meshkit-docs/examples.html
[online]: http://sigma.mcs.anl.gov/meshkit-library/
[CGM build instructions]: http://sigma.mcs.anl.gov/cgm/building-cgm/
[MOAB build instructions]: http://sigma.mcs.anl.gov/moab/building-moab/
[Lasso build instructions]: http://sigma.mcs.anl.gov/lasso/building-lasso/