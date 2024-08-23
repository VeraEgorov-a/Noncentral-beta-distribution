% Example: Computing the Noncentral Beta Distribution using Bpqxy
clear vars

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


