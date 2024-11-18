# Two Compartment Model for Brain Clearance

This repository contains code for simulating a three-compartment model for brain clearance. In this model

Compartment 1 represents Brain (B)

Compartment 2 represents CSF (C)

Compartment 3 represents Plasma (P)

## Model Description

The model consists of the following parameters:

- A: Amyloid production during wake
- sigma_A: Scaling factor to denote a loss in productioin during sleep
- $\text r_{bc}$: Amyloid transfer from Brain to CSF during wake
- $\text sigma_{bc}$: Scaling factor to denote increase in amyloid transfer from Brain to CSF during sleep
- $\text r_{bp}: Amyloid transfer from Brain to Plasma during wake
- $\text sigma_{bp}$: Scaling factor to denote increase in amyloid transfer from Brain to Plasma during sleep
- $\text r_{cp}$: Amyloid transfer from CSF to Plasma during wake
- $\text sigma_{cp}$: Scaling factor to denote increase in amyloid transfer from CSF to Plasma during sleep
- $\text r_{p}$: Amyloid transfer from Plasma to Systemic Circulation during wake
- $\text sigma_{p}$: Scaling factor to denote increase in amyloid transfer from Plasma to Systemic Circulation during sleep

## Optimised Parameters Values

| Parameter | Value |
|-----------|-------|
| A       | 12   |
| $\text sigma_{A}$  | 0.8   |
| $\text r_{bc}$ | 1.5   |
| $\text sigma_{bc}$   | 2.5    |
| $\text r_{bp}$   | $\text{r}_{bp} = \frac{\text{r}_{bc} \cdot (1 - 133 \cdot \text{r}_{cp})}{133 \cdot \text{r}_{cp}}$|
| $\text sigma_{bp}         | 3.99     |
| $\text r_{cp}         | 3.99     |
| $\text sigma_{cp}         | 0.005     |
| $\text r_{p}         | 0.28     |
| $\text sigma_{p}         | 2.58     |

## Steady State Solutions

The steady state solutions for the model with the `a21` parameter are as follows:

- B_w(t) = A_wake / (k - a21)
- C_w(t) = (k * A_wake) / (a12_wake * (k - a21))
- B_s(t) = A_sleep / (k - a21)
- C_s(t) = (k * A_sleep) / (a12_sleep * (k - a21))

## Usage

To use the code, follow the instructions in the respective script files.

