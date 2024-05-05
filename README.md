# Two Compartment Model for Brain Clearance

This repository contains code for simulating a two-compartment model for brain clearance. In this model

Compartment 1 represents CSF (C)

Compartment 2 represents Blood (B)

## Model Description

The model consists of the following parameters:

- a21: Amyloid transfer from Blood to CSF
- a12_wake: Amyloid transfer from CSF to blood during wake state
- a12_sleep: Amyloid transfer from CSF to blood during sleep state
- A_wake: Amyloid production during wake
- A_sleep: Amyloid production during sleep
- k: Amyloid clearance from blood

## Parameters Values

| Parameter | Value |
|-----------|-------|
| a21       | 0.3   |
| a12_wake  | 0.1   |
| a12_sleep | 0.5   |
| A_wake    | 69    |
| A_sleep   | 1     |
| k         | 2     |

## Steady State Solutions

The steady state solutions for the model with the `a21` parameter are as follows:

- B_w(t) = A_wake / (k - a21)
- C_w(t) = (k * A_wake) / (a12_wake * (k - a21))
- B_s(t) = A_sleep / (k - a21)
- C_s(t) = (k * A_sleep) / (a12_sleep * (k - a21))

## Usage

To use the code, follow the instructions in the respective script files.

