# Numerical Algorithms

This repository contains implementations of various numerical algorithms used for solving mathematical problems. Each algorithm is implemented in C++ and includes examples for easy understanding and usage.

## Implemented Algorithms

- **Approximation (`Approximation.cpp`)**  
	DefinItoin of a class named ChebyshevApproximation that does an approximation of a function using Chebyshev polinomials, as well as derivation and integration based on this approximation.
  

- **Extrapolation (`Extrapolator.cpp`)**  
  Definition of a global function Limit for numerical finding of limit values using Richardson extrapolation.
  

- **Local Minimum Finder (`FindLocalMinimum.cpp`)**  
  Algorithms to numerically calculate local minimum of one-variable functions.
  

- **Numerical Integration (`Integrations.cpp`)**  
  Implementation of three global functions that support the numerical integration RombergIntegration, TanhSinhIntegration, AdaptiveIntegration
  

- **Interpolation (`Interpolators.cpp`)**  
	Definition of six classes that cover numerical interpolation: LinearInterpolator, PolynomialInterpolator, SplineInterpolator (cubic spline), TrigonometricInterpolator (Gauss), BarycentricInterpolator (barycentric rational interpolation). All of these classes are derived from an abstract base class AbstractInterpolator to allow polymorphism
	
- **Lossy Compression and Decompression (`LossyCompressDecompress.cpp`)**  
  Functions written to support lossy compression of slowly varying sequences, based on the discrete Fourier transformation, and for approximate reconstruction of the original sequence from the compressed sequence.

- **Matrix Operations (`MatrixOperations.cpp`, `MatrixOperationsV2.cpp`)**  
  Functions for performing various operations on matrices.

- **Runge-Kutta 4th Order Integrator (`RK4Integrator.cpp`)**  
  Implementation of the 4th order Runge-Kutta method for solving regular first-order differential equations using function RK4Integrator, and solving a system of first-order differential equations using the function RK4SystemIntegrator.

- **Root-Finding (`ZeroOfFunctions.cpp`)**  
  A collection of functions to support numerical calculation of nonlinear equations: BracketRoot - enclose a root of the function using the "bracketing" technique, RegulaFalsiSolve - user cna choose type of algorithm: RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic}, RiddersSolve, NewtonRaphsonSolve: for this algorithm the root doesn't have to be bracketed, PolyRoots: for polynomial functions.

