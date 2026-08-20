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

# Bounds for the reciprocal

"""
Note: From this point onward, we abuse notation slightly and write t for log(t).
"""

## Generalised framework

def a(r,t):
    return 1 + fdiv(log(1 + r*exp(-t)),t)
def b(r,t):
    return 1 + fdiv(log(a(r,t)),log(t))
def A(omega,t,xi1,xi2):
    if omega > 0 and t >= exp(1):
        return 70.6995
    if omega == 0 and t >= log(3):
        return 58.096*sqrt(1 + 9*exp(-2*t))*(1 + fdiv(log(1 + 9*exp(-2*t)),2*t))
def B(omega,t,xi1,xi2):
    return fdiv(2,3) + 4.43795*power(omega,fdiv(3,2))*power(log(t),fdiv(1 - 3*xi2,2))*power(t,1-fdiv(3*xi1,2))
def A3(d,t,xi1,xi2):
    return fdiv(power(A(0,t,xi1,xi2)*(1 + fdiv(log(2),t)),fdiv(1,6)),power(d,0.75))*exp(fdiv(3*g0()*d,4*power(t,fdiv(2,3))*power(log(t),fdiv(1,3))))
def A4(d,omega,t,r,xi1,xi2):
    pwr = 4.43795*power(omega,fdiv(3,2))*power(fdiv(power(log(t),1 - xi2),power(t,xi1)),fdiv(3,2))
    term1 = power(1 + r*exp(-t),pwr)
    return term1*A3(d,t,xi1,xi2)*A(omega,t-r,xi1,xi2)*power(a(r,t),B(omega,t,xi1,xi2))
def ub(t,d,omega,W,xi1,xi2,mode=1):
    # Upper bounds for the logarithic derivatives..
    R = fdiv(d + omega*log(t),power(t,xi1)*power(log(t),xi2))
    # First, check whether the condition holds..
    if 1 < fdiv(fdiv(omega,d)*log(t) + 1,2*(power(a(R,t),-xi1)*power(b(R,t),-xi2) + 1)):
        check = True
    else:
        check = False
    if mode == 1:
        beta = fdiv(fdiv(1,d*W) + 1,power(a(R,t),-xi1)*power(b(R,t),-xi2) + 1)
        lamda = fdiv(8*beta,1 - beta)*(1 + fdiv(d,omega*log(t)))*power(1 - fdiv(d + 2*d*power(a(R,t),-xi1)*power(b(R,t),-xi2),omega*log(t)),-1)*power(omega - fdiv(d + 2*d*power(a(R,t),-xi1)*power(b(R,t),-xi2),log(t)),-1)
        o1 = lamda*(B(omega,t,xi1,xi2) + fdiv(3*xi1,4) + fdiv(1,6) + fdiv(3*xi2*log(log(t)),4*log(t)) + fdiv(log(A4(d,omega,t,R,xi1,xi2)),log(t))) + fdiv(1+beta,d*(1-beta)) 
        o2 = fdiv(1,fdiv(1,W) + 2*d)
        return max(o1,o2), round(beta,10), check
    if mode == 2:
        beta = fdiv(1,power(a(R,t),-xi1)*power(b(R,t),-xi2) + 1)
        lamda = fdiv(8*beta,1 - beta)*(1 + fdiv(d,omega*log(t)))*power(1 - fdiv(d + 2*d*power(a(R,t),-xi1)*power(b(R,t),-xi2),omega*log(t)),-1)*power(omega - fdiv(d + 2*d*power(a(R,t),-xi1)*power(b(R,t),-xi2),log(t)),-1)
        o1 = lamda*(B(omega,t,xi1,xi2) + fdiv(3*xi1,4) + fdiv(1,6) + fdiv(3*xi2*log(log(t)),4*log(t)) + fdiv(log(A4(d,omega,t,R,xi1,xi2)),log(t))) + fdiv(1+beta,d*(1-beta)) 
        o2 = fdiv(1,2*d)
        return max(o1,o2), round(beta,10), check
def pairings(t,W=54,J=200,increments=0.001,split_at=2,xi1=1,xi2=0):
    pairs = []
    if xi1 == 1 and xi2 == 0:
        Z = 5.558691
        om = 3.4812610455
    elif xi1 == 1 and xi2 == -1:
        Z = 21.233
        om = 1.0144550824
    else:
        Z = 53.989
        om = 1.5046945258
    for V in [W + i*increments for i in range(int(split_at/increments))]:
        a, b, c = ub(t,fdiv(1,Z),om,V,xi1,xi2)
        pairs.append((V,a))
    for V in [W + 0.5*i for i in range(int(2*split_at),int(2*J))]:
        a, b, c = ub(t,fdiv(1,Z),om,V,xi1,xi2)
        pairs.append((V,a))
    return pairs
def R(t,W=54,J=200,xi1=1,xi2=0,incs=0.0001):
    if xi1 == 1 and xi2 == 0:
        con = 24.303
        d = 0.00049
    elif xi1 == 1 and xi2 == -1:
        con = 110.6
        d = 0.00069
    else:
        con = 216.667
        d = 0.00346104
    VQ_pairs = pairings(t,W,J,incs,2,xi1,xi2)
    Y = fdiv(VQ_pairs[-1][1],VQ_pairs[-1][0])
    j = 1
    while j < J:
        if fdiv(1,VQ_pairs[j-1][0]) - fdiv(1,VQ_pairs[j][0]) < 0:
            print(j, VQ_pairs[j-1][0], VQ_pairs[j-1][0], fdiv(1,VQ_pairs[j-1][0]) - fdiv(1,VQ_pairs[j][0]))
        Y += VQ_pairs[j-1][1]*(fdiv(1,VQ_pairs[j-1][0]) - fdiv(1,VQ_pairs[j][0]))
        j += 1
    """
    If required, we choose d optimially to 8dp as follows:
    d = 0.000000001 
    for v in [0.0001,0.00001,0.00000001]:
        while A3(d+v,t,xi1,xi2)*max(1,exp(Y + 216.667*(d+v))) < A3(d,t,xi1,xi2)*max(1,exp(Y + con*d)):
            d += v
        d -= v
    d += v
    """
    return A3(d,t,xi1,xi2)*max(1,exp(Y + con*d)), d, VQ_pairs[0][1]

## Classical bounds 

for W in [10,8,7,6,5.56]:
    if W > 5.56:
        m1, m2, m3 = R(log(H()),W,500,1,0,0.00001)
        print(f" ${W}$ & ${print_pretty(int(roundup(m1,0)),5,False)}$ & ${round(m2,8)}$ \\\\")
    else:
        m1, m2, m3 = R(log(H()),5.56,500,1,0,0.0000001)
        print(f" ${W}$ & ${print_pretty(int(roundup(m1,0)),5,False)}$ & ${round(m2,8)}$ \\\\")

"""
Output:
 $10$ & $678$ & $0.00049$ \\
 $8$ & $686$ & $0.00049$ \\
 $7$ & $701$ & $0.00049$ \\
 $6$ & $793$ & $0.00049$ \\
 $5.56$ & $1237$ & $0.00049$ \\
"""

## Littlewood bounds

for W in [30,25,22,21.5,21.24]:
    if W > 21.24:
        m1, m2, m3 = R(log(H()),W,500,1,-1,0.00001)
        print(f" ${W}$ & ${print_pretty(int(roundup(m1,0)),5,False)}$ & ${round(m2,8)}$ \\\\")
    else:
        m1, m2, m3 = R(log(H()),W,500,1,-1,0.000001)
        print(f" ${W}$ & ${print_pretty(int(roundup(m1,0)),5,False)}$ & ${round(m2,8)}$ \\\\")

"""
Output:
 $30$ & $721$ & $0.00069$ \\
 $25$ & $731$ & $0.00069$ \\
 $22$ & $790$ & $0.00069$ \\
 $21.5$ & $942$ & $0.00069$ \\
 $21.24$ & $1980$ & $0.00069$ \\
"""

## Korobov--Vinogradov bounds (\sigma < 1)

for W in [70,60,55,54,53.99]:
    if W > 53.99:
        m1, m2, m3 = R(log(H()),W,500,fdiv(2,3),fdiv(1,3),0.00001)
        print(f" ${W}$ & ${print_pretty(int(roundup(m1,0)),5,False)}$ & ${round(m2,8)}$ \\\\")
    else:
        m1, m2, m3 = R(log(H()),W,500,fdiv(2,3),fdiv(1,3),0.000001)
        print(f" ${W}$ & ${print_pretty(int(roundup(m1,0)),5,False)}$ & ${round(m2,8)}$ \\\\")

"""
Output:
 $70$ & $467$ & $0.00346104$ \\
 $60$ & $473$ & $0.00346104$ \\
 $55$ & $488$ & $0.00346104$ \\
 $54$ & $5223$ & $0.00346104$ \\
 $53.99$ & $6368$ & $0.00346104$ \\
"""

## Korobov--Vinogradov bounds (\sigma \geq 1)

def optimised_ub_sg1(t,Z=5.558691,xi1=1,xi2=0):
    omga = 0.5
    for v in [10**(-p) for p in range(11)]:
        while ub(t,fdiv(1,Z),omga + v,W,xi1,xi2,2)[0] < ub(t,fdiv(1,Z),omga,W,xi1,xi2,2)[0]:
            omga += v
        omga -= v
    omga += v
    comp1, comp2, check = ub(t,fdiv(1,Z),omga,W,xi1,xi2,2)
    return comp1
def proc_sg1(Z=5.558691,xi1=1,xi2=0):
    T0 = 5
    for v in [1000,100,10,1,0.1,0.01]:
        while 24.303*power(fdiv(T0+v,log(T0+v)),fdiv(1,3)) < optimised_ub_sg1(T0+v,Z,xi1,xi2):
            T0 += v
    return T0, optimised_ub_sg1(T0,Z,xi1,xi2)

proc_sg1(53.989,fdiv(2,3),fdiv(1,3))

"""
Output:
(6186.0700000000015,
 mpf('216.666947007816352835602272862041947555235890331876968509374450767324729790079419447795142710742787437134600499418127680453'))
"""
