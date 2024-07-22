M = Manifold(4, 'M', structure='Lorentzian')
X.<x0,r,th,ph> = M.chart('x0 r:(0,oo) th:(0,pi):\\theta ph:period=2*pi:\\phi')
a = function('a')
k, c = var('k c')

g = M.metric()
g[0,0] = -1
g[1,1] = a(x0/c)^2 / (1 - k*r^2)
g[2,2] = a(x0/c)^2 * r^2
g[3,3] = a(x0/c)^2 * r^2 * sin(th)^2

Ric = g.ricci()
R = g.ricci_scalar()
LHS = Ric - 1/2 * R * g
LHS.set_name('LHS', '\\mathrm{LHS}')

ddx0 = X.frame()[0]
rho = var('rho', latex_name='\\rho')
q, G = var('q, G')
Lam = var('Lam', latex_name='\\Lambda')

# Stress-energy tensor for perfect fluid
T_upp = (rho*c^2 + q) * ddx0 * ddx0 + q * g.inverse()
T_mix = T_upp.down(g, 1)
T_low = T_mix.down(g, 0)

RHS = 8*pi*G/c^4 * T_low - Lam * g
RHS.set_name('RHS', '\\mathrm{RHS}')

LHS_trace = LHS.up(g, 0).trace(0, 1)
RHS_trace = RHS.up(g, 0).trace(0, 1)

# Spatial slice at time t
t = var('t')
N = Manifold(3, 'N', ambient=M, structure="Riemannian", start_index=1)
Y = N.chart('r:(0,oo) th:(0,pi):\\theta ph:period=2*pi:\\phi')
embed = N.diff_map(M, {(Y, X): [c*t, r, th, ph]})
proj_N = M.diff_map(N, {(X, Y): [r, th, ph]})
proj_t = M.scalar_field({X: x0/c})
N.set_embedding(embed, inverse=proj_N, var=t, t_inverse={t: proj_t})
