# Task list

## General

- [] Add support to package `sf` and drop `sp` out of the `Imports` list.
- [] Try adding support to package `data.table` as a means of improving memory usage and computation speed.

## scheduleSPSANN

- [] `initial.acceptance`. Pass a vector with two numeric values defining the minimum and maximum initial acceptance probability, i.e. `initial.acceptance = c(0.95, 0.99)`. This is needed so that users won't pass an unnecessarily large value to `initial.temperature`.
- [] `initial.temperature`. Ideally, `spsann` would estimate the initial temperature internally so that users do not have to worry about finding an appropriate value to reach the desired `initial.acceptance`.
