/*!
 \page algorithm_assygen Assembly Geometry and Mesh Script Generator

<b>Name:</b> %AssyGen

<b>External dependencies:</b> CGM --with a geometry engine

<b>Input:</b> %AssyGen input file

<b>Output:</b> Assembly geometry file (.sat or .stp) and Cubit meshing scripts for meshing this assembly

<b>Interface(s) used:</b> iGeom

<b>Setup:</b>

Need to call PrepareIO with with command line arguments to input the AssyGen input file. 

<b>Notes:</b>

Default test file case: test/algs/test_assygen.cpp uses input file: data/assygen_default.inp. To run other AssyGen input files using the program use executable 'test_assygen' in 'MeshKit/test/algs' folder; specying the filename:- 'Your Dir> test_assygen s1' where, s1.inp is the AssyGen input file.
See paper:
http://link.springer.com/chapter/10.1007/978-3-642-15414-0_21
*/
