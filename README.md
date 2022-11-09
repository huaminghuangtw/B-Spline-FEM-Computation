B-Spline-FEM-Computation
========================

> FE code that discretizes a two-dimensional Poisson’s equation using B-Spline basis functions.

### Setup
1. Recursively clone this project to your machine: `git clone --recursive https://github.com/huaminghuangtw/B-Spline-FEM-Computation.git`  
2. Use Cmake to build/generate project files corresponding to your preferred complier/IDE, e.g., Visual Studio 2019  
   (Remember to tick the `PYBIND11_INSTALL` box, otherwise there will be `ImportError: DLL not found` message)
   <p align="left">
      <img width="500" src="figures/cmake.png"/>
   </p>
3. Open project solution file (in this case `splinecomputation.sln`)  
4. Go to Solution Explorer > Right-click on `INSTALL` > Select `build`
   <p align="left">
      <img width="280" src="figures/install.png"/>
   </p>
5. A new folder called `install` will be generated under your build folder  
6. After this step you can go to the `install` folder and run the _splinekernel testrunner_ or _any python scripts_
   <p align="left">
      <img width="600" src="figures/installfolder.png"/>
   </p>

---

### How to use the code?

In the **install** folder, you could play around with:

- splinekernel testrunner
  - CMD  
    <img align="left" width=600 src="figures/cmd_test.png"><br clear="all"/>
  - Git Bash  
    <img align="left" width="600" src="figures/gitbash_test.png"><br clear="all"/>
  - Visual Studio
    * Go to Solution Explorer > Right-click on `splinekernel_testrunner` > Select `Debug` > Select `Start New Instance`
    * Go to Solution Explorer > Right-click on `splinekernel_testrunner` > Select `Set as StartUp Project` > Click `F5`
- Python scripts
  ```python
  python plotBSplineBasis.py
  ```
  <img align="left" width="400" src="figures/plotBSplineBasis.py.png"><br clear="all"/>
  
  ```python
  python plotBSplineCurve.py
  ```
  <img align="left" width="400" src="figures/plotBSplineCurve.py.png"><br clear="all"/>
  
  ```python
  python plotBSplineBasis2D.py
  ```
  <img align="left" width="400" src="figures/plotBSplineBasis2D.py.png"><br clear="all"/>
  
  ```python
  python plotBSplineSurface.py
  ```
  <img align="left" width="400" src="figures/plotBSplineSurface.py.png"><br clear="all"/>
  
  ```python
  python laplaceProblem.py
  ```
  <img align="left" width="400" src="figures/laplaceProblem.py.png"><br clear="all" />

---

### Contact
If you have any question or suggestion, feel free to contact me at huaming.huang.tw@gmail.com. Contributions are also welcomed. Please open a [pull-request](https://github.com/huaminghuangtw/B-Spline-FEM-Computation/compare) or an [issue](https://github.com/huaminghuangtw/B-Spline-FEM-Computation/issues/new) in this repository.
