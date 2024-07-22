var('a')
dS = Manifold(4, 'dS', structure='Lorentzian')
X.<t,r,th,ph> = dS.chart(f't r:(0,{a}) th:(0,pi) ph:period=2*pi')

g = dS.metric()
g[0,0] = -(1 - r^2/a^2)
g[1,1] = 1/(1 - r^2/a^2)
g[2,2] = r^2
g[3,3] = r^2*sin(th)^2

ddt,ddr,ddth,ddph = X.frame()[:]
nabla = g.connection()
Ric = g.ricci()
R = g.ricci_scalar()