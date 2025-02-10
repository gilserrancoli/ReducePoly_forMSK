# Reduction of muscle-tendon kinematic parameterization
This repository contains the main function to reduce the polynomials used to parameterize muscle-tendon lengths and moment arms.
The method presented here is the one published in:
Harba, M.; Badia, J.; Serrancolí, G. XXXX. YYYY.


In this code, Main_Polynomials.m first compute full polynomials parameterizing muscle-tendon lengths (lMT) and moment arms from joint coordinates, following the method presented in XXX and also used in other optimal control problems YYY. First, random positions of the skeleton are created covering the whole range of motion that we expect to have within the movements. 

The added value presented in the mentioned publication is the reduction of those polynomials (within Polynomial_reduction.m, which is called at the end of Main_Polynomials.m). First, the polynomial coefficients are ordered according to the significance of this term within the whole polynomial. To do this, the standard error (SE) was computed, the t-term was computed dividing the coefficient by SE, and then the p-value associated with it was obtained. All coefficients were sorted by its computed p-value. Then, starting from the first coefficient, we test if the resulting polynomial meets the requiring accuracy criteria (measured by the maximum of RMSE of the lMT and moment arms). If the resulting polynomial meets the criteria, that polynomial would be the "reduced polynomial", if not, we add the following coefficient of the list and we asses the accuracy again. See the following flowchart.

![following flowchart](img/flowchart.png).
