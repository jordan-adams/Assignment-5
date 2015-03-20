# Homework 5 Problem 3
# Jordan Adams

from math import exp
from numpy import pi,cos,tan,linspace,ones,copy

# Imports for required math functions

def gauss1(n):
    a = linspace(3,4*n-1,n)/(4*n+2)
    
    x = cos(pi*a+1/(8*(n**2)*tan(a)))

# User defined gaussian function (1)
    
    epsilon = 1e-15
    
    delta = 1
    
    while delta>epsilon:

# This while block handles the operation for the desired values        
        
        p1 = ones(n,float)
        p2 = copy(x)
        for k in range(1,n):
            p1,p2 = p2,((2*k+1)*x*p2-k*p1)/(k+1)
        dp = (n+1)*(p1-x*p2)/(1-(x**2))
        dx = p2/dp
        x -= dx
        delta = max(abs(dx))

    w = 2*(n+1)*(n+1)/((n**2)*(1-x*x)*(dp**2))
    
    # This calculates the weight values

    return x,w

def gauss2(n,a,b):
    x,w = gauss1(n)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

# User defined gaussian function (2)

def f(c): 
    if c<1e-4:
        return c**2+7*c**3/2+91*c**4/12+13*c**5+13859*c**6/720+2309*c**7/90
    else:
        return c**3/((1-c)**5*(exp(c/(1-c))-1))
h=30

# Number of slices

a=0
b=1

# Limits of integration

x,w=gauss2(h,a,b)
s=0
for k in range(h):
    s+=w[k]*f(x[k])

print("Value of integral is {0:.2f}".format(s))

# Outputs the integral value for reading
# The program asked to compare it against a known value, but to solve 
# W=(boltzmanconstant)(T**4) I would need a temperature to solve for it9io