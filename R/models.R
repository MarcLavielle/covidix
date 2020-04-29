sir1 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, L0, tmax}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
L_0 = L0

tm = min(t, tmax)
ks = k0-a1*tm
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_L = kd*I - L/d
ddt_D = L/d
")

sir2 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, L0, tmax, tau1, a2}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
L_0 = L0

tm = min(t, tmax)
if tm < tau1
   ks = k0 - a1*tm
else
   ks = k0 - a1*tau1 + a2*(tm-tau1)
end
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_L = kd*I - L/d
ddt_D = L/d
")

sir3 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, L0, tmax, tau1, tau2, a0, a2}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
L_0 = L0

tm = min(t, tmax)
if tm < tau1
   ks = k0 - a0*(tau1-tm)
elseif tm < tau2
   ks = k0 - a1*(tm-tau1)
else
   ks = k0 - a1*(tau2-tau1) + a2*(tm-tau2)
end
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_L = kd*I - L/d
ddt_D = L/d
")


