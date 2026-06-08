---
title: Modeling World Cup Traffic Flow
subtitle: Numerical Computation Final Project Spring 26
author:
- Rose Lancaster
- JP Appel
date: June 9, 2026
---

# Intro

## Quantum Graph

* associate a length and (partial) differential equation to each edge

## Continuity/Advection Equation

$$\delta_t \rho + \nabla \cdot (\textbf{u}\rho) = \sigma$$

* $\rho$ is distribution of quantity
* $\textbf{u}$ is vector field describing velocity
* $\textbf{u}\rho$ is called flux in many physical settings (electromagnetism)
* $\sigma$ is a source or sink function for the quantity in concern over time

# Scenarios

## Advection Dirichlet B.C.

### Equation Description

$$\delta_t P + \delta_x (c P) = 0$$

::: notes
* use method of lines
* solve via finite differences
* also solvable by a different technique (finite volumes)
* if c is independent of spatial variable (x) then this becomes a problem we have already worked with
:::

### Initial Conditions

![**PLACEHOLDER FOR INITIAL CONDITIONS**](pictures/advection_initcond.png)

### Velocity function

![**PLACEHOLDER FOR PLOT OF C**](pictures/advection_velocity.png)

:::notes
* here the velocity function is a flipped hat
* given this, we expect a wave to slow down in the region less than 1
:::

### Animation

![**PLACEHOLDER FOR ANIMATION**](pictures/test.png)

::: notes
* this animation shows what we expect, wave moves along with a near constant speed until it hits the middle third
:::

## Advection Periodic B.C.

### Periodic Boundary Conditions

:::notes
* we can introduce periodic boundary conditions to make the problem behave as a closed system
:::

![**PLACEHOLDER FOR ANIMATION**](pictures/test.png)

## Two node Conservative

![**PLACEHOLDER FOR ANIMATION**](pictures/test.png)

## Two node Non-conservative

![**PLACEHOLDER FOR ANIMATION**](pictures/test.png)

## Multi node

![**PLACEHOLDER FOR ANIMATION**](pictures/test.png)

# Modeling

## Public Dataset
