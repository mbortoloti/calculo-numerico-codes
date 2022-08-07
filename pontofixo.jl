#
#   Este código aproxima a raíz de uma equação do tipo f(x) = 0.
#
using Printf;

# Code parameters
ϵ = 1.e-6;
maxiter = 1000;

# Function setting

# Example 01
    f(x) = log(x) - x + 2;

    # Iteration functions
    #ϕ₁(x) = log(x) + 2;
    ϕ₂(x) = exp(x-2);

#Example 02
 #   f(x) = cos(x) - 3x;

    # Iteration functions
#    ϕ₁(x) = cos(x) / 3.0;

# Initial guess
x = 0.4;
iter = 0;

while true
    global x
    global iter
    afx = abs(f(x));
    @printf("%5d %12.8e\n",iter,afx);
    
    if afx < ϵ
        @printf("Solution was achieved!");
        break;
    end

    iter = iter + 1;
    if iter > maxiter
        @printf("Maximum of iterations was achieved! Stopping...");
        break;
    end

    # Create new sequence term
    x = ϕ₂(x);
end
