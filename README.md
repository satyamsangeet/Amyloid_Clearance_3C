# Three Compartment Model for Brain Clearance

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
- $\text r_{bp}$: Amyloid transfer from Brain to Plasma during wake
- $\text sigma_{bp}$: Scaling factor to denote increase in amyloid transfer from Brain to Plasma during sleep
- $\text r_{cp}$: Amyloid transfer from CSF to Plasma during wake
- $\text sigma_{cp}$: Scaling factor to denote increase in amyloid transfer from CSF to Plasma during sleep
- $\text r_{p}$: Amyloid transfer from Plasma to Systemic Circulation during wake
- $\text sigma_{p}$: Scaling factor to denote increase in amyloid transfer from Plasma to Systemic Circulation during sleep

## Optimised Parameters Values

| Parameter | Value |
|-----------|-------|
| A       | 12.063   |
| $\text sigma_{A}$  | 0.782   |
| $\text r_{bc}$ | 1.623   |
| $\text sigma_{bc}$   | 2.505    |
| $\text r_{bp}$   | 0.199    |
| $\text sigma_{bp}$         |  4.618    |
| $\text r_{cp}$         | 0.00572     |
| $\text sigma_{cp}$         | 5.190     |
| $\text r_{p}$         | 0.300     |
| $\text sigma_{p}$         | 4.670     |

## Usage

To use the code, follow the instructions in the respective script files.

