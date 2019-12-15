# README #

### Summary ###
Following the paper:

Gao Y, Bouix S. Statistical Shape Analysis using 3D Poisson
Equation---A Quantitatively Validated Approach. Medical Image
Analysis. 2016 Jan 15.

this repository contains code to perform statistical analysis on 3D
shapes with arbitrary topology. It contains an iterative 3D Poisson
equation solver using ITK (VNL) sparse system.

### What is this repository for? ###

* Statistical shape analysis

### How do I get set up? ###

* CMake (cmake.org)
* ITK (itk.org, Compilation needed)
* alglib library
* GNU gsl

### Docker

The built version is also available on Docker: https://hub.docker.com/r/gaoyi/spom/

```bash
docker pull gaoyi/spom
```

### How to use:
Give two sets of shapes, this software/algorithm performs statistal analysis on them. It output the average/mean shape of the input shapes. Moreover, on the mean shape, the locations where the two groups have significant difference are highlighted.

* Compile the code so you have the executable file named "mainSumOfTwoPoissonShapeAnalysis"
* Save the two groups of shapes to be analyzed as binary volumetric images. Could be in the format of nifti, nrrd, mhd/mha, etc.
* Put the full path/names of the image names of one shape group into a text file, say list1.txt
* Put the full path/names of the image names of the other shape group into a text file, say list2.txt
* Run 
  * If the input shapes have ALREADY been registered, run:
```bash
./path-to/mainSumOfTwoPoissonShapeAnalysis list1.txt list2.txt meanImageName.nrrd meanVTPName.vtp 0
```
  * If the input shapes have NOT been registered, 
```bash
./path-to/mainSumOfTwoPoissonShapeAnalysis list1.txt list2.txt meanImageName.nrrd meanVTPName.vtp 1
```

### Who do I talk to? ###

* Yi Gao (gaoyi@gatech.edu)

### Usage agreement ###

In scientific publications (journal publications, conference papers,
technical reports, presentations at conferences and meetings) that use
this code, you should cite the following paper:

Gao Y, Bouix S. Statistical Shape Analysis using 3D Poisson
Equation---A Quantitatively Validated Approach. Medical Image
Analysis. 2016 Jan 15.

### To do ###

