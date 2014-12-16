function g = LieEl(phi)

T = Lg();
Tphi = T*phi;
g = expm(Tphi);
