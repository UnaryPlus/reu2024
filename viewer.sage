from scipy.integrate import solve_ivp
from PIL import Image
import line_profiler

def debug(msg, a, b):
  print(msg + ':', a)
  return b

def schwarzschild_metric(rs):
  M = Manifold(4, 'M', structure='Lorentzian')
  X.<t,r,th,ph> = M.chart(f't r:({rs},oo) th:(0,pi) ph:period=2*pi')

  g = M.metric()
  g[0,0] = -(1 - rs/r)
  g[1,1] = 1/(1 - rs/r)
  g[2,2] = r^2
  g[3,3] = r^2*sin(th)^2
  
  return M

def sph_to_cart(r, th, ph):
  return r * sin(th) * cos(ph), r * sin(th) * sin(ph), r * cos(th)

def sphere(r_c, th_c, ph_c, radius):
  r_c, th_c, ph_c, radius = float(r_c), float(th_c), float(ph_c), float(radius)
  x_c, y_c, z_c = sph_to_cart(r_c, th_c, ph_c)
  def termination_condition(t, r, th, ph):
    x, y, z = sph_to_cart(r, th, ph)
    return (x - x_c)^2 + (y - y_c)^2 + (z - z_c)^2 - radius^2
  return termination_condition

def geodesic_integrator(M, termination_conditions):
  dim = M.dim()
  v_test = M.tangent_space(M.an_element()).an_element()
  c_test = M.integrated_geodesic(M.metric(), (SR.var('s'), 0, 1), v_test)
  
  X = M.default_chart()
  coords = list(X)
  velocities = X.symbolic_velocities()
  des = velocities + c_test.system()[0]
  des = [ fast_callable(de, vars=(coords + velocities), domain=float) for de in des ]
  eqn = lambda s, y: [ de(*y) for de in des ]

  # List comprehension does not work here (in a hard to catch way). I HATE PYTHON
  # Also why does mapping over a list not return a list?
  events = list(map(lambda f: lambda s, y: f(*y[:dim]), termination_conditions)) 
  for e in events: e.terminal = True
  
  def integrate(v, s_max):
    p = v.parent().base_point()
    init = [ float(x) for x in X(p) + tuple(v) ]
    return solve_ivp(eqn, (0, s_max), init, events=events) # 78% of total time

  return integrate

M = schwarzschild_metric(1)
g = M.metric()
ddt, ddr, ddth, ddph = M.default_frame()[:]

p_obs = M((0, 3, pi/2, 0))
dot = g.at(p_obs)
basis = [ ddt, -ddr, ddph, ddth ] # time, forward, right, down
basis = [ v.at(p_obs) for v in basis ]
basis = [ v / sqrt(abs(dot(v, v))) for v in basis ]
B_obs = M.tangent_space(p_obs).basis('B', from_family=basis)

spheres = [ sphere(3, pi/2, 5*pi/4, 1), sphere(0, 0, 0, 1.0001) ]
integrate = geodesic_integrator(M, spheres)

def local_color(v, s_max):
  incoming_light = v - sqrt(dot(v,v)) * basis[0]
  result = integrate(incoming_light, s_max)
  if result.status == 0: 
    # s_max was reached
    return b'\x00'
  elif result.status == 1:
    if len(result.t_events[-1]) > 0: 
      # geodesic hit black hole
      return b'\x11'
    else:
      # geodesic hit a light source
      return b'\xff'
  else:
    raise Exception('An error occurred during integration.')

WIDTH = 50

def pixel_to_vector(i, j):
  # return value is not necessarily a unit vector
  # division by WIDTH is to make integration step size reasonable
  z, x, y = WIDTH/4, i - WIDTH/2 + 0.5, j - WIDTH/2 + 0.5
  return (z * basis[1] + x * basis[2] + y * basis[3]) / WIDTH

@line_profiler.profile
def main():
  colors = []
  s_max = 40
  for row in range(WIDTH):
    for col in range(WIDTH):
      v = pixel_to_vector(col, row)
      colors.append(local_color(v, s_max))
  img = Image.frombytes('L', (WIDTH, WIDTH), b''.join(colors))
  img.save('black_hole.png')

main()