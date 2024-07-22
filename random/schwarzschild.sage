var('t r th ph x y z rs')

Sch = Manifold(4, 'Sch', structure='Lorentzian')
X = Sch.chart('t r:(0,oo) th:(0,pi) ph:period=2*pi')

g = Sch.metric()
g[0,0] = -(1 - rs/r)
g[1,1] = 1/(1 - rs/r)
g[2,2] = r^2
g[3,3] = r^2*sin(th)^2

ddt,ddr,ddth,ddph = X.frame()[:]
nabla = g.connection()

Y = Sch.chart('t x y z', coord_restrictions=(x != 0, y != 0, z != 0))
sph_to_cart = X.transition_map(Y, [t, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
sph_to_cart.set_inverse(t, sqrt(x^2 + y^2 + z^2), atan2(sqrt(x^2 + y^2), z), atan2(y, x))

SS = Manifold(3, 'SS', ambient=Sch, structure="Riemannian", start_index=1)
X_ = SS.chart('r:(0,oo) th:(0,pi) ph:period=2*pi')
embed = SS.diff_map(Sch, {(X_, X): [0, r, th, ph]})
proj = Sch.diff_map(SS, {(X, X_): [r, th, ph]})
SS.set_embedding(embed, inverse=proj)
