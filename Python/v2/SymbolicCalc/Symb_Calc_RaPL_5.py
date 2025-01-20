from sympy import symbols
from sympy import solve
from sympy import diff
from sympy import exp

def diff_conv(exp,var):
    return str(diff(exp,var))

# Define variables
th4 = symbols("th[3]") 
th5 = symbols("th[4]") 
th6 = symbols("th[5]") 
th7 = symbols("th[6]") 
s = symbols("s")

t = symbols("t")
t0 = symbols("tb[0]")
t1 = symbols("tb[1]")
t2 = symbols("tb[2]")
t3 = symbols("tb[3]")
t4 = symbols("tb[4]")
t5 = symbols("tb[5]")

A = symbols("A")
V = symbols("V")
al = symbols("al")

S = ((t1-t0)*exp(th4) + \
     (t2-t1)*(exp(th4)+exp(th5)) + \
     (t3-t2)*(exp(th5)+exp(th6)) + \
     (t4-t3)*(exp(th6)+s) + \
     (t5-t4)*(s+exp(th7)))/2 + \
      exp(th7) / al

   
txtf = open("SymbolicCalc\Results_RaPL_5.txt",mode="w")

# Solve for s
print("Solution for s",file=txtf)
sol = solve(S - A, s, dict=True)
print("s =",sol[0][s],file=txtf)

# Calculate derivatives
print("\nDerivatives",file=txtf)
s = sol[0][s]
Ra = exp(th4)/(t1-t0)*(t-t0) / V
print("if tb[0] <= t <= tb[1]:",file=txtf)
print("\tdFdTh[0, nf] =",diff(Ra,th4),file=txtf)
Ra = (exp(th4) + (exp(th5)-exp(th4))/(t2-t1)*(t-t1)) / V
print("elif tb[1] <= t <= tb[2]:",file=txtf)
print("\tdFdTh[0, nf] =",diff(Ra,th4),file=txtf)
print("\tdFdTh[0, nf+1] =",diff(Ra,th5),file=txtf)
Ra = (exp(th5) + (exp(th6)-exp(th5))/(t3-t2)*(t-t2)) / V
print("elif tb[2] <= t <= tb[3]:",file=txtf)
print("\tdFdTh[0, nf+1] =",diff(Ra,th5),file=txtf)
print("\tdFdTh[0, nf+2] =",diff(Ra,th6),file=txtf)
Ra = (exp(th6) + (s-exp(th6))/(t4-t3)*(t-t3)) / V
print("elif tb[3] <= t <= tb[4]:",file=txtf)
print("\tdFdTh[0, nf] =",diff(Ra,th4),file=txtf)
print("\tdFdTh[0, nf+1] =",diff(Ra,th5),file=txtf)
print("\tdFdTh[0, nf+2] =",diff(Ra,th6),file=txtf)
print("\tdFdTh[0, nf+3] =",diff(Ra,th7),file=txtf)
Ra = (s + (exp(th7)-s)/(t5-t4)*(t-t4)) / V
print("elif tb[4] <= t <= tb[5]:",file=txtf)
print("\tdFdTh[0, nf] =",diff(Ra,th4),file=txtf)
print("\tdFdTh[0, nf+1] =",diff(Ra,th5),file=txtf)
print("\tdFdTh[0, nf+2] =",diff(Ra,th6),file=txtf)
print("\tdFdTh[0, nf+3] =",diff(Ra,th7),file=txtf)
Ra = exp(th7) * exp(-al*(t-t5)) / V
print("else:",file=txtf)
print("\tdFdTh[0, nf+3] =",diff(Ra,th7),file=txtf)

txtf.close()