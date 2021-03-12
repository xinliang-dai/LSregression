# LSregression

### Introduction

- Implementing _Gauss-Newton_ with different _conjugate gradient_ methods, including

  - without _conjugate gradient_

  - _standard conjugate gradient_

  - _preconditioned conjugate gradient_

  - _CG-Steihaug_

  - _preconditioned CG-Steihaug_

- Implementing _Levenberg-Marquardt_ with some methods to increase speed and robustness

  - damping strategy of lambda, due to Nielsen:

    - [Methods for nonlinear least squares problems](http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf) - introduction and comparison

    - [Damping Parameter in Marquardt's Method](http://www.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf) - original & statistic detail

  - _CG-Steihaug_

  - _preconditioning_
 

### Getting Start

Carrying out data-fitting test: [nonlinear_datafitting_test.m](../main/nonlinear_datafitting_test.m)

Nonlinear Least-Squares solver: [pcg_steihaug_gauss_newton.m](../main/nls_solvers/pcg_steihaug_gauss_newton.m)

### ClassDef

Example cases saved in [data folder](../main/data)

 - [dim_1_case](../main/data/dim_1_case.m):
f_i(x) = exp(d_i*x)

 - [dim_2_case.m](../main/data/dim_2_case.m):
f_i(x) = exp( d_i * x_1 ) * sin( d_i * x_2 )

 - [dim_2_quadratic_case.m](../main/data/dim_2_quadratic_case.m):
 f_i(x) = d_i * x_1  +  d_i.^2 * x_2
