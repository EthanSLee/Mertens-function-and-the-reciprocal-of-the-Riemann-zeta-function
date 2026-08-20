import mpmath as mpm
from mpmath import mpf, sqrt, pi, log, exp, power, fdiv, ceil, floor, euler
from numpy import linspace

mpm.mp.dps = 120

# Preamble

def print_pretty(value,dp=2,override=False):
    if override == False:
        exponent = 7
    else:
        exponent = 4
    if value < 10**exponent:
        if isinstance(value, int) == True:
            return f"{value}"
        else:
            return f"{roundup(value,dp)}"
    else:
        while value / (10 ** exponent) >= 10:
            exponent += 1
        return f"{roundup(value / (10 ** exponent),dp)} \cdot 10^{{{exponent}}}"
def roundup(y,m):
    return round(ceil(y*power(10,m))*power(10,-m),m)
def g0():
    return euler
def H():
    return 3000175332800

# Bounds for Mertens function

"""
Note: From this point onward, we abuse notation slightly and write t for log(t).
"""

## Generalised framework

def F(x,W,l1,l2):
    return power(fdiv(x,W),fdiv(1,l1+1))*power(fdiv(l1+1,log(x) - log(W)),fdiv(l2,l1+1))
def condition1(x,W,l1,l2):
    if l2 >= 0:
        LHS = F(x,W,l1,l2)
    else:
        LHS = power(1 - fdiv(log(W) + l2*log(fdiv(l1+1,log(x) - log(W))),log(x)),-fdiv(l2,l1+1))*F(x,W,l1,l2)
    if LHS > log(H()):
        return True
    else:
        return False
def condition2(x,W,l1,l2,l3,l4,R1,R2):
    LHS = fdiv(R1*power(log(H()),l3 + 1)*power(log(log(H())),l4),pi*(l3 + 1))
    if l2 >= 0:
        RHS = 2*exp(1) + fdiv(R1,pi)*(fdiv(exp(1)*power(log(H()),l3 - l1 - 1)*power(log(log(H())),l4 - l2),W) + fdiv(power(log(H()),l3)*power(log(log(H())),l4),H()*x)) + exp(power(fdiv(x,W),fdiv(1,l1+1)) - fdiv(x,4))*(fdiv(R2*power(log(H()),2),2*pi) + fdiv(power(4,2)*56,3*pi) + 4*exp(1 + fdiv(1,x))*exp(-fdiv(3*x,4)))
    else:
        RHS = 2*exp(1) + fdiv(R1,pi)*(fdiv(exp(1)*power(log(H()),l3 - l1 - 1)*power(log(log(H())),l4 - l2),W) + fdiv(power(log(H()),l3)*power(log(log(H())),l4),H()*x)) + exp(power(fdiv(x*log(x),W),fdiv(1,l1+1)) - fdiv(x,4))*(fdiv(R2*power(log(H()),2),2*pi) + fdiv(power(4,2)*56,3*pi) + 4*exp(1 + fdiv(1,x))*exp(-fdiv(3*x,4)))
    if LHS > RHS:
        return True
    else:
        return False

## Classical bounds

def R_ves(W):
    if W == 10:
        return 678
    if W == 8:
        return 686
    if W == 7:
        return 701
    if W == 6:
        return 793
    if W == 5.56:
        return 1237
def c1(x,W,R1):
    return 2*exp(1) + fdiv(12*R1,23*pi*W)*power(fdiv(x,W),-fdiv(1,24)) + fdiv(4*exp(2 + fdiv(1,x)),sqrt(W*x))

for k in [5.56,6,7,8,10]:
    if condition2(k*power(log(H()),2),k,1,0,fdiv(11,12),0,R_ves(k),exp(146)) == False:
        print(f"Disaster at W = {k}")
    else:
        print(f"${k}$ & ${roundup(k*power(log(H()),2),2)}$ & ${R_ves(k)}$ & ${roundup(c1(k*power(log(H()),2),k,R_ves(k)),2)}$ \\\\")

"""
Output:
$5.56$ & $4589.2$ & $1237$ & $33.56$ \\
$6$ & $4952.38$ & $793$ & $22.21$ \\
$7$ & $5777.77$ & $701$ & $18.16$ \\
$8$ & $6603.17$ & $686$ & $16.34$ \\
$10$ & $8253.96$ & $678$ & $14.06$ \\
"""

## Littlewood bounds

def R_ves(W):
    if W == 30:
        return 721
    if W == 25:
        return 731
    if W == 22:
        return 790
    if W == 21.5:
        return 942
    if W == 21.24:
        return 1980
def c2(x,W,R1):
    return 2*exp(1) + fdiv(12*R1,23*pi*W)*power(log(H()),-fdiv(1,12))*power(log(log(H())),fdiv(1,4)) + fdiv(4*exp(2 + fdiv(1,x)),sqrt(fdiv(W*x,log(x))))

for k in [21.24,21.5,22,25,30]:
    j = 100
    for increment in [10**(2-a) for a in range(5)]:
        while condition1(j+increment,k,1,-1) == False:
            j += increment
    if condition2(j,k,1,-1,fdiv(11,12),-fdiv(3,4),R_ves(k),exp(146)) == False:
        print(f"Disaster at W = {k}")
    else:
        print(f"${k}$ & ${round(j,2)}$ & ${R_ves(k)}$ & ${roundup(c2(j,k,R_ves(k)),2)}$ \\\\")

"""
Output:
$21.24$ & $10352.03$ & $1980$ & $21.48$ \\
$21.5$ & $10490.06$ & $942$ & $13.08$ \\
$22$ & $10755.87$ & $790$ & $11.73$ \\
$25$ & $12360.29$ & $731$ & $10.57$ \\
$30$ & $15066.77$ & $721$ & $9.66$ \\
"""

## Korobov--Vinogradov bounds

def R_ves(W):
    if W == 70:
        return 467
    if W == 60:
        return 473
    if W == 55:
        return 488
    if W == 54:
        return 5223
    if W == 53.99:
        return 6368
def c3(x,W,R1):
    return 2*exp(1) + fdiv(3*R1,5*pi*W)*power(log(log(H())),-fdiv(1,12)) + fdiv(4*exp(2 + fdiv(1,x)),power(W,fdiv(3,5))*power(x,fdiv(2,5)))

for k in [53.99,54,55,60,70]:
    j = 5
    for increment in [10**(2-a) for a in range(5)]:
        while condition1(exp(j+increment),k,fdiv(2,3),fdiv(1,3)) == False:
            j += increment
    if condition2(exp(j),k,fdiv(2,3),fdiv(1,3),fdiv(2,3),fdiv(1,4),R_ves(k),exp(146)) == False:
        print(f"Disaster at W = {k}")
    else:
        print(f"${k}$ & ${round(j,2)}$ & ${R_ves(k)}$ & ${roundup(c3(exp(j),k,R_ves(k)),2)}$ \\\\")

"""
Output:
$53.99$ & $10.01$ & $6368$ & $25.85$ \\
$54$ & $10.01$ & $5223$ & $22.19$ \\
$55$ & $10.03$ & $488$ & $7.02$ \\
$60$ & $10.11$ & $473$ & $6.85$ \\
$70$ & $10.27$ & $467$ & $6.63$ \\
"""

## Threshold computations

def bound1(x):
    return fdiv(2.9189,power(x,2))
def bound2(x):
    return 33.56*exp(-sqrt(fdiv(x,5.56)))
def bound3(x):
    return 21.48*x*exp(-sqrt(fdiv(x*(log(x) - log(21.24)),2*21.24)))
def bound4(x):
    return 25.85*x*exp(-power(fdiv(x,53.99),fdiv(3,5))*power(fdiv(5,3*(log(x) - log(53.99))),fdiv(1,5)))

yA = 100.0
for increment in [10**(2-a) for a in range(5)]:
    while bound1(yA + j) < bound2(yA + j):
        yA += j
print(yA)
yB = yA
for increment in [10**(2-a) for a in range(5)]:
    while bound2(yB + j) < bound3(yB + j):
        yB += j
print(yB)
yC = 5.0
for increment in [10**(1-a) for a in range(4)]:
    while bound3(exp(yC + j)) < bound4(exp(yC + j)):
        yC += j
print(yC)

"""
Output:
1650.7699999999977
138529.3300000094
25.539999999999996
"""
