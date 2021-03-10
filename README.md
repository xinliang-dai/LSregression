# LSregression

### Introduction

Implementing _Gauss-Newton_ with different _conjugate gradient_ methods, including

- without _conjugate gradient_

- _standard conjugate gradient_

- _preconditioned conjugate gradient_

- _CG-Steihaug_

- _preconditioned CG-Steihaug_


### Getting Start

Carrying out data-fitting test: [nonlinear_datafitting_test.m](../tree/main/nonlinear_datafitting_test.m)

Nonlinear Least-Squares solver: [pcg_steihaug_gauss_newton.m](../tree/main/nls_solvers/pcg_steihaug_gauss_newton.m)

### ClassDef

Example cases saved in [data folder](../tree/main/data)

 - [dim_1_case](../tree/main/data/dim_1_case.m):
f_i(x) = exp(d_i*x)

 - [dim_2_case.m](../tree/main/data/dim_2_case.m):
f_i(x) = exp( d_i * x_1 ) * sin( d_i * x_2 )

 - [dim_2_quadratic_case.m](../tree/main/data/dim_2_quadratic_case.m):
 f_i(x) = d_i * x_1  +  d_i.^2 * x_2
