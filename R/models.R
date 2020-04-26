sir1 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
;ks = k0*exp(-a1*t)
ks = k0-a1*t
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_D = kd*delay(I,d)
")

sir2 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, tau1, a2}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
if t < tau1
  ks = k0 - a1*t
else
  ks = k0 - a1*tau1 + a2*(t-tau1)
end
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_D = kd*delay(I,d)
")

sir3 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, tau1, tau2, a0, a2}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
if t < tau1
  ks = k0 - a0*(tau1-t)
elseif t < tau2
  ks = k0 - a1*(t-tau1)
else
  ks = k0 - a1*(tau2-tau1) + a2*(t-tau2)
end
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_D = kd*delay(I,d)
")


