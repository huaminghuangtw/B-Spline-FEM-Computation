B-Spline-FEM-Computation
=======================

> FE code that discretizes a two-dimensional Poisson’s equation using B-Spline basis functions.

---

## Setup
1. Recursively clone this project to your machine: `git clone --recursive git@github.com:hmhuang0501/B-Spline-FEM-Computation.git`
2. Use Cmake to build/generate project files corresponding to your preferred complier/IDE, e.g., Visual Studio 2019  
   (Remember to tick the "PYBIND11_INSTALL" box, otherwise there will be "ImportError: DLL not found" message)  
   ![cmake](figures/cmake.png)
3. Open project solution file (in this case **splinecomputation.sln**)
4. Go to Solution Explorer, right-click on "INSTALL", choose "build"  
   ![install](figures/install.png)
5. A new folder called "install" will be generated under your build folder
6. After this step you can go to the "install" folder and run the _splinekernel testrunner_ or _any python scripts_  
   ![installfolder](figures/installfolder.png)

---

## How To Use

In the **install** folder, you can play around with:

- splinekernel testrunner
  - CMD  
    ![cmd_test](figures/cmd_test.png)
  - Git Bash  
    ![gitbash_test](figures/gitbash_test.png)
- python scripts
  ```python
  python plotBSplineBasis.py
  ``` 
  ![plotBSplineBasis.py](figures/plotBSplineBasis.py.png)  
  
  ```python
  python plotBSplineCurve.py
  ```
  ![plotBSplineCurve.py](figures/plotBSplineCurve.py.png)  
  
  ```python
  python plotBSplineBasis2D.py
  ```
  ![plotBSplineBasis2D.py](figures/plotBSplineBasis2D.py.png)  
  
  ```python
  python plotBSplineSurface.py
  ```
  ![plotBSplineSurface.py](figures/plotBSplineSurface.py.png)  
  
  ```python
  python laplaceProblem.py
  ```
  ![laplaceProblem.py](figures/laplaceProblem.py.png)  

---
