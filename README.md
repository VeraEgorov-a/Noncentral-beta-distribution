# Bpqxy

`Bpqxy` is a MATLAB function designed to compute the Noncentral Beta Distribution. 

## Function Signature

```matlab
[Bpq, ierr] = Bpqxy(x, y, p, q)
```

## Inputs

- **`x`**: Noncentrality parameter (lambda, should be positive)
- **`y`**: The quantile or value at which the noncentral beta CDF is evaluated (should be in the range [0, 1])
- **`p`**: Shape parameter \( p \) (must be positive)
- **`q`**: Shape parameter \( q \) (must be positive)

## Outputs

- **`Bpq`**: Computed value of the noncentral beta distribution
- **`ierr`**: Error flag indicating the success or issues in computation:
  - `0`: Computation successful
  - `1`: Some loss of accuracy is expected in the computed function value
  - `2`: Computation failed due to one or more input values being out of bounds

## Description

The `Bpqxy` function computes the noncentral beta distribution $B_{p,q}(x, y)$ using various numerical techniques. The function is designed to handle a wide range of input values and provides an error flag to indicate if there were any issues during the computation.

The noncentral beta distribution function, denoted as $B_{p,q}(x, y)$, is defined by the following series expansion:

$$
B_{p,q}(x, y) = e^{-\frac{x}{2}} \sum_{j=0}^{\infty} \frac{1}{j!} \left(\frac{x}{2}\right)^j I_{y}(p + j, q),
$$

where $I_{y}(p, q)$ is the central beta distribution function (implemented in the `betaincreg(y,p,q)` function), $p$ and $q$ are the shape parameters, and $x$ is the noncentrality parameter. 


## Usage

Here is a basic example of how to use the `Bpqxy` function:

```matlab
% Example: Computing the Noncentral Beta Distribution using Bpqxy

% Given parameters
x = 4.5;  % Noncentrality parameter (lambda)
y = 0.12;  % Input value (should be in the range [0, 1])
p = 3.0;  % Shape parameter p (must be positive)
q = 4.0;  % Shape parameter q (must be positive)

% Call the Bpqxy function
[I,ierr] = Bpqxy(x, y, p, q);

% Call the MATLAB built-in function
nu1=2*p;
nu2=2*q;
w=nu2*y/(nu1*(1-y));
Im=ncfcdf(w,nu1,nu2,x);

% Display the results
fprintf('Computed noncentral beta distribution value (Bpqxy): I = %.8e\n', I);
fprintf('Computed noncentral beta distribution value (MATLAB): I = %.8e\n', Im);
```

## Authors

V. Egorova, A. Gil, J. Segura and N. M. Temme


## Contact

For questions or comments, please contact Amparo Gil at [amparo.gil@unican.es].

