sir0 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, L0}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
W_0 = W0
D_0 = D0
L_0 = L0
ks = k0*exp(-a1*t)
ddt_I = (ks - kr - kd)*I
ddt_W = ks*I
ddt_L = kd*I -L/d
ddt_D = L/d
")

sir1 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, d, kd, kr, I0, W0, D0, L0, AW, AD, phi}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
Wa_0 = W0
Da_0 = D0
L_0 = L0
ks = k0*exp(-a1*t)
ddt_I = (ks - kr - kd)*I
ddt_Wa = ks*I
ddt_L = kd*I -L/d
ddt_Da = L/d

;W = Wa*(1 + AW*cos(6.28*t/7-phi))
;D = Da*(1 + AD*cos(6.28*t/7-phi))
W = Wa + AW*cos(6.28*t/7-phi)*I*ks
D = Da + AD*cos(6.28*t/7-phi)*L/d
")


sir2 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, a2, d, kd, kr, I0, W0, D0, L0, AW, AD, phi, tau1}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
Wa_0 = W0
Da_0 = D0
L_0 = L0
if t< tau1
   ks = k0 + a2*t
else
   A = (a2*tau1 + k0)/exp(-a1*tau1)
   ks = A*exp(-a1*t)
end
ddt_I = (ks-kr)*I - kd*I
ddt_Wa = ks*I
ddt_L = kd*I -L/d
ddt_Da = L/d

;W = Wa*(1 + AW*cos(6.28*t/7-phi))
;D = Da*(1 + AD*cos(6.28*t/7-phi))
W = Wa + AW*cos(6.28*t/7-phi)*I*ks
D = Da + AD*cos(6.28*t/7-phi)*L/d
")

sir3 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, a2, a3, d, kd, kr, I0, W0, D0, L0, AW, AD, phi, tau1, tau2}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
Wa_0 = W0
Da_0 = D0
L_0 = L0
if t<tau1
   ks=k0 + a2*t
elseif t<tau2
   ks=k0 + a2*t + a3*(t-tau1)
else
   A = (a2*tau2 + a3*(tau2-tau1) + k0)/exp(-a1*tau2)
   ks = A*exp(-a1*t)
end
ddt_I = (ks-kr)*I - kd*I
ddt_Wa = ks*I
ddt_L = kd*I -L/d
ddt_Da = L/d

;W = Wa*(1 + AW*cos(6.28*t/7-phi))
;D = Da*(1 + AD*cos(6.28*t/7-phi))
W = Wa + AW*cos(6.28*t/7-phi)*I*ks
D = Da + AD*cos(6.28*t/7-phi)*L/d
")

sir4 <- inlineModel("
[LONGITUDINAL]
input = {k0, a1, a2, a3, a4, d, kd, kr, I0, W0, D0, L0, AW, AD, phi, tau1, tau2, tau3}
EQUATION:
odeType=stiff
t0 = 0
I_0 = I0
Wa_0 = W0
Da_0 = D0
L_0 = L0
if t<tau1
   ks=k0 + a2*t
elseif t<tau2
   ks=k0 + a2*t + a3*(t-tau1)
elseif t<tau3
   ks=k0 + a2*t + a3*(t-tau1) + a4*(t-tau2)
else
   A = (a2*tau3 + a3*(tau3-tau1) + a4*(tau3-tau2) + k0)/exp(-a1*tau3)
   ks = A*exp(-a1*t)
end
ddt_I = (ks-kr)*I - kd*I
ddt_Wa = ks*I
ddt_L = kd*I -L/d
ddt_Da = L/d

;W = Wa*(1 + AW*cos(6.28*t/7-phi))
;D = Da*(1 + AD*cos(6.28*t/7-phi))
W = Wa + AW*cos(6.28*t/7-phi)*I*ks
D = Da + AD*cos(6.28*t/7-phi)*L/d
")
