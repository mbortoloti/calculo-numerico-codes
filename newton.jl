#
#   Newton code for solving f(x) = 0
#

using Printf;

# Tolerance setting
ϵ = 1.e-12;

maxiter = 1000;


# function setting
#f(x) = 2.0 * (x-1.0)^3 * cos(x) * exp(x);
#f(x) = (x-1)^7;
f(x) = cos(pi * x / 2.0)*exp(-1.0/x^2);

# Derivative of for
#df(x) = 2*(x-1.0)^2*exp(x)*(x*cos(x)-x*sin(x)+2*cos(x)+sin(x));
#df(x) = 7*(x-1)^6;
df(x) = -exp(-1/x^2)*(pi*sin(pi*x/2)*x^3-4*cos(pi*x/2))/(2*x^3);

# Iteration function
ϕ(x) = x - f(x) / df(x);

#Exact Solution
xstar = 1.0; 

# Initial guess
xk = 0.1;

# Iteration counter
k = 0;

# Only for convergence order
ekp1 = NaN;
ek = NaN;
ekm1 = NaN;

while true

    global xk,k,ek,ekm1,ekp1;

    nfx = abs(f(xk));

    #Error calculation
    ekm1 = ek;
    ek = ekp1;
    ekp1 = abs(xk - xstar);
    α = log(ekp1/ek) / log(ek/ekm1);
    # print info
    @printf("%5d    %10.8e    %10.8e    %10.8e    %10.8f\n",k,xk,nfx,ek,α);


    if nfx < ϵ
        println("Solution was found!");
        break;
    end

    # Update iteration counter
    k = k + 1;

    if k > maxiter
        println("Maximum of iterations was achieved! Stopping...");
        break;
    end

    # Calculate new sequence term
    xk = ϕ(xk);

end


