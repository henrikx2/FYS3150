---
title: "Numerical integration methods"
date: "21.10.2019"
author: "Henrik Lind Petlund | *github.com/henrikx2*"
---

## Abstract
Key results, why is this work worthwhile? Credibility for my claims.

## 1 Introduction
The ground state corrolation energy between two electrons in a helium atom can be determined by solving a 6-dimensional integral. This integral is derived by modelling the wave function of each electron as an single-particle wave function of the electron in the hydrogen atom. For an electron *i* in the 1s state, the dimensionless and unnormalized single-particle wave function can be expressed as

\[ \psi_{1s}({\bf r}_i) = e^{-\alpha r_{i}} \]

where $\alpha$ is a parameter, and

\[ {\bf r}_i =  x_i {\bf e}_x + y_i {\bf e}_y +z_i {\bf e}_z \]


\[ {r}_i = \sqrt{{x}_i^2+{y_i^2}+z_i^2} \]

The parameter $\alpha=2$ gives the charge of the helium atom ($Z=2$). Further, the wave function of the two 1s electrons are given by

\[ \Psi(r_{1}+r_{2}) = e^{-\alpha (r_{2}+r_{2})} \]

The integral to solve in this report is the expectation value for the corrolation energy between the two electrons in the helium atom. The corrolation energy depends on the classical Coluomb interactions of the two electrons, and is given by

\[\begin{align} \langle \frac{1}{|{\bf r}_1-{\bf r_2}|} \rangle = \int_{-\infty}^{\infty}d{\bf r}_1d{\bf r}_2e^{-2\alpha({\bf r}_1+{\bf r}_2)}\frac{1}{{|\bf r}_1-{\bf r}_2|} \end{align} \label{equ} \tag{1}\]

This (unnormalized) integral can be solved on closed form to be $5\pi^2/16^2\approx0.19276571$. This can be showed to be correct by splitting the integral into $r_1$ and $r_2$ dependance, and finding these integrals in an encyclopedia (i.e. Rottman).

## 2 Methods

### 2.1 Gauss Quadrature
Gauss Quadrature is a method that uses orthogonal polynomials to estimate the integral. This methdos are referenced in [@gaussQuad] and is quite extensively to cover for this report, and is therefore just sited here.

#### 2.2.1 Gauss-Legendre Quadrature


#### 2.2.2 Gauss-Laguerre Quadrature (Improved Gauss Quadrature)

### 2.2 Monte Carlo Integration

#### 2.2.1 Brute force Monte Carlo Integration

#### 2.2.2 Improved Monte Carlo Integration

#### 2.2.3 Improved Monte Carlo Integration with Parallization

## 3 Resulsts

## 4 Discusson

## 5 Conclusion

## 6 Appendix

## Rererences
