# Reduction of muscle-tendon kinematic parameterization
This folder contains examples of tracking simulations to evaluate the difference of using full polynomials to parameterize muscle-tendon lengths and moment arms and the reduced polynomials, presented in:
Harba, M.; Badia, J.; Serrancolí, G. XXXX. YYYY.

The main function is <a href="OCP_GC\RunAllSimulations.m">RunAllSimulations.m (from OCP_GC folder)</a>. By default, the script runs tracking optimal-control simulations of experimental gait trials using a musculoskeletal model with a smooth right-knee contact model representing a knee prosthesis. It evaluates four gait conditions—normal, metatarsophalangeal (MTP), treadmill, and bouncy—and, for each condition, solves the same optimal control problem under two surrogate parameterizations of musculotendon geometry: (1) full polynomials and (2) reduced polynomials for muscle–tendon lengths and moment arms. This yields eight simulations in total.
