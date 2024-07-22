M = Manifold(4, 'M', structure='Lorentzian')
X.<t,r,th,ph> = M.chart('t r:(0,oo) th:(0,pi) ph:period=2*pi')
a = function('a')
k = var('k')

g = M.metric()
g[0,0] = -1
g[1,1] = a(t)^2 / (1 - k*r^2)
g[2,2] = a(t)^2 * r^2
g[3,3] = a(t)^2 * r^2 * sin(th)^2

Ric = g.ricci()
R = g.ricci_scalar()
G = Ric - 1/2 * R * g
G.set_name('G')
