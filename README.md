# easyPlace
Reimplementation of the electrostatic-based VLSI placement algorithm: ePlace and ePlace-MS with clean C++ code. This code is for those who try to understand the electrostatics-based placement. 

Reimplementation of Abacus is applied for internal std cell legalization.

Independent set matching, local reordering and global swap are applied for internal detailed placement.

# How to Build
Requires CMake 3.16+ and a C++17 compiler. macOS Apple Clang and Linux GCC are supported. OpenMP is enabled when available; otherwise the program runs without OpenMP. Placement images are saved as BMP files without requiring X11.

From the repository root:
```sh
cmake -S . -B build
cmake --build build --target ePlace -j 8
```

# How to Run
Run the included synthetic example (100 standard cells, 4 fixed pins and 184 nets):
```sh
./build/main/ePlace -aux ./examples/smoke/smoke.aux -targetDensity 0.8 -targetOverflow 0.1 -internalLegal 1 -internalDP -outputPath ./output
```
This runs quadratic placement, global placement, internal Abacus legalization and detailed placement. Results are saved under `output/smoke/`, including `smoke-eDP.pl` and `Graphs/Detailed placement result.bmp`. This example verifies the standard-cell flow; it is not a benchmark of large mixed-size designs.

For an existing Bookshelf benchmark, keep its `.nodes`, `.nets`, `.wts`, `.pl` and `.scl` files alongside the `.aux` file, then run:
```sh
./build/main/ePlace -aux /path/to/adaptec4.aux -targetDensity 1.0 -fullPlot -targetOverflow 0.1 -legalizerPath /path/to/legalizer -outputPath ./output
```
The external legalizer directory must contain a compatible `ntuplace3` executable. For standard-cell designs, `-internalLegal 1 -internalDP` uses the internal implementation instead. The program recreates the output subdirectory named after the benchmark on each run.

# Options
* -aux: specify input aux file
* -targetDensity: specify target density
* -targetOverflow: specify target overflow
* -legalizerPath: specify the path of external legalizer (support ntuplace3 only for now) and call ntuplace3 to complete legalization and detailed placement
* -outputPath: specify the output path 
* -fullPlot: plot cell locations during the placement process
* -noQP: skip initial placement(quadratic placement)
* -nomGP: skip mixed-size global placement
* -nomLG: skip macro legalization and cell global placement
* -nocGP: skip cell global placement
* -noLegal: skip legalization
* -addNoise: add random noise before mGP

# References
### Global Placement:
* [ePlace](https://dl.acm.org/doi/10.1145/2699873)
* [ePlace-MS](https://ieeexplore.ieee.org/document/7008518)
* [RePlAce](https://ieeexplore.ieee.org/document/8418790)
* [Xplace](https://github.com/cuhk-eda/Xplace)

### Legalization:
* [Abacus](https://dl.acm.org/doi/10.1145/1353629.1353640)

### Detailed placement techniques:
* [NTUplace3](https://ieeexplore.ieee.org/document/4544855)
* [ABCDPlace](https://ieeexplore.ieee.org/document/8982049)
* [FastPlace](https://ieeexplore.ieee.org/document/1560039)

# Authors
Ziang Ge and Yikai Liu, supervised by Prof. [Pingqiang Zhou](https://faculty.sist.shanghaitech.edu.cn/faculty/zhoupq/home.html) at [ShanghaiTech University](https://www.shanghaitech.edu.cn/)
