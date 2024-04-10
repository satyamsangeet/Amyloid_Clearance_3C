# Two Compartment Model for Brain Clearance

This repository contains code for simulating a two-compartment model for brain clearance. In this model

Compartment 1 represents Blood (B)

Compartment 2 represents CSF (C)

## Model Description

The model consists of the following parameters:

- a12: Amyloid transfer from Blood to CSF
- a21_wake: Amyloid transfer from CSF to blood during wake state
- a21_sleep: Amyloid transfer from CSF to blood during sleep state
- A_wake: Amyloid production during wake
- A_sleep: Amyloid production during sleep
- k: Amyloid clearance from blood

## Parameters Values

| Parameter | Value |
|-----------|-------|
| a12       | 0.3   |
| a21_wake  | 0.1   |
| a21_sleep | 0.5   |
| A_wake    | 69    |
| A_sleep   | 1     |
| k         | 2     |

## Steady State Solutions

The steady state solutions for the model with the `a12` parameter are as follows:

- B_w(t) = A_wake / k
- C_w(t) = ((a12 + k) * A_wake) / (a21_wake * k)
- B_s(t) = A_sleep / k
- C_s(t) = ((a12 + k) * A_sleep) / (a21_sleep * k)

The steady state solutions for the model without the `a12` parameter are as follows:

- B_w(t) = A_wake / k
- C_w(t) = A_wake / a21_wake
- B_s(t) = A_sleep / k
- C_s(t) = A_sleep / a21_sleep

## Usage

To use the code, follow the instructions in the respective script files.

