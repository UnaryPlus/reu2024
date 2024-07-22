

# This file was *autogenerated* from the file viewer.sage
from sage.all_cmdline import *   # import sage library

_sage_const_4 = Integer(4); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3); _sage_const_1p0001 = RealNumber('1.0001'); _sage_const_50 = Integer(50); _sage_const_0p5 = RealNumber('0.5'); _sage_const_40 = Integer(40)
from scipy.integrate import solve_ivp
from PIL import Image

def debug(msg, a, b):
  print(msg + ':', a)
  return b

def schwarzschild_metric(rs):
  M = Manifold(_sage_const_4 , 'M', structure='Lorentzian')
  X = M.chart(f't r:({rs},oo) th:(0,pi) ph:period=2*pi', names=('t', 'r', 'th', 'ph',)); (t, r, th, ph,) = X._first_ngens(4)

  g = M.metric()
  g[_sage_const_0 ,_sage_const_0 ] = -(_sage_const_1  - rs/r)
  g[_sage_const_1 ,_sage_const_1 ] = _sage_const_1 /(_sage_const_1  - rs/r)
  g[_sage_const_2 ,_sage_const_2 ] = r**_sage_const_2 
  g[_sage_const_3 ,_sage_const_3 ] = r**_sage_const_2 *sin(th)**_sage_const_2 
  
  return M

def sph_to_cart(r, th, ph):
  return r * sin(th) * cos(ph), r * sin(th) * sin(ph), r * cos(th)

def sphere(r_c, th_c, ph_c, radius):
  r_c, th_c, ph_c, radius = float(r_c), float(th_c), float(ph_c), float(radius)
  x_c, y_c, z_c = sph_to_cart(r_c, th_c, ph_c)
  def termination_condition(t, r, th, ph):
    x, y, z = sph_to_cart(r, th, ph)
    return (x - x_c)**_sage_const_2  + (y - y_c)**_sage_const_2  + (z - z_c)**_sage_const_2  - radius**_sage_const_2 
  return termination_condition

def geodesic_integrator(M, termination_conditions):
  dim = M.dim()
  v_test = M.tangent_space(M.an_element()).an_element()
  c_test = M.integrated_geodesic(M.metric(), (SR.var('s'), _sage_const_0 , _sage_const_1 ), v_test)
  
  X = M.default_chart()
  coords = list(X)
  velocities = X.symbolic_velocities()
  des = velocities + c_test.system()[_sage_const_0 ]
  des = [ fast_callable(de, vars=(coords + velocities), domain=float) for de in des ]
  eqn = lambda s, y: [ de(*y) for de in des ]

  # List comprehension does not work here (in a hard to catch way). I HATE PYTHON
  # Also why does mapping over a list not return a list?
  events = list(map(lambda f: lambda s, y: f(*y[:dim]), termination_conditions)) 
  for e in events: e.terminal = True
  
  def integrate(v, s_max):
    p = v.parent().base_point()
    init = [ float(x) for x in X(p) + tuple(v) ]
    return solve_ivp(eqn, (_sage_const_0 , s_max), init, events=events)

  return integrate

M = schwarzschild_metric(_sage_const_1 )
g = M.metric()
ddt, ddr, ddth, ddph = M.default_frame()[:]

p_obs = M((_sage_const_0 , _sage_const_3 , pi/_sage_const_2 , _sage_const_0 ))
dot = g.at(p_obs)
basis = [ ddt, -ddr, ddph, ddth ] # time, forward, right, down
basis = [ v.at(p_obs) for v in basis ]
basis = [ v / sqrt(abs(dot(v, v))) for v in basis ]
B_obs = M.tangent_space(p_obs).basis('B', from_family=basis)

spheres = [ sphere(_sage_const_3 , pi/_sage_const_4 , pi, _sage_const_1 ), sphere(_sage_const_0 , _sage_const_0 , _sage_const_0 , _sage_const_1p0001 ) ]
integrate = geodesic_integrator(M, spheres)

def local_color(v, s_max):
  incoming_light = v - sqrt(dot(v,v)) * basis[_sage_const_0 ]
  result = integrate(incoming_light, s_max)
  if result.status == _sage_const_0 : 
    # s_max was reached
    return b'\x00'
  elif result.status == _sage_const_1 :
    if len(result.t_events[-_sage_const_1 ]) > _sage_const_0 : 
      # geodesic hit black hole
      return b'\x11'
    else:
      # geodesic hit a light source
      return b'\xff'
  else:
    raise Exception('An error occurred during integration.')

WIDTH = _sage_const_50 

def pixel_to_vector(i, j):
  # return value is not necessarily a unit vector
  # division by WIDTH is to make integration step size reasonable
  z, x, y = WIDTH/_sage_const_4 , i - WIDTH/_sage_const_2  + _sage_const_0p5 , j - WIDTH/_sage_const_2  + _sage_const_0p5 
  return (z * basis[_sage_const_1 ] + x * basis[_sage_const_2 ] + y * basis[_sage_const_3 ]) / WIDTH

colors = []
s_max = _sage_const_40 
for row in range(WIDTH):
  for col in range(WIDTH):
    v = pixel_to_vector(col, row)
    colors.append(local_color(v, s_max))
img = Image.frombytes('L', (WIDTH, WIDTH), b''.join(colors))
img.save('black_hole.png')
