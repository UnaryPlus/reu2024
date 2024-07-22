from sage.manifolds.differentiable.pseudo_riemannian import PseudoRiemannianManifold

def default(x, y):
    if x == None: return y
    else: return x

def arc_length(c, g, low=None, high=None):
    """
    Length of curve `c` (from `low` to `high`) under Riemannian metric `g`.
    By default, `low` and `high` are set to the endpoints of the curve.
    """
    I = c.domain()
    t = I.canonical_coordinate()
    low = default(low, I.inf())
    high = default(high, I.sup())

    v = c.tangent_vector_field()
    speed = sqrt(g.along(c)(v, v))
    return speed.expr().integral(t, low, high)

def einstein(g, name='G'):
    """
    Einstein tensor of pseudo-Riemannian metric `g`.
    """
    G = g.ricci() - 1/2 * g.ricci_scalar() * g
    G.set_name(name)
    return G

class HyperbolicSpace(PseudoRiemannianManifold):
    """
    `n`-dimensional hyperbolic space equipped with the standard half-space and unit-ball charts.
    (The unit-ball chart is not created until `ball_chart()` is called.)
    By default, `name` is set to `'H^n'`.
    """
    def __init__(self, n, name=None, latex_name=None, metric_name=None, metric_latex_name=None, start_index=1, metric_constant=1):
        name = default(name, f'H^{n}')
        latex_name = default(latex_name, r'\mathbf{H}^{' + str(n) + '}')
        PseudoRiemannianManifold.__init__(self, n, name, latex_name=latex_name, metric_name=metric_name, metric_latex_name=metric_latex_name, signature=n, start_index=start_index)
        
        xs = [f'x{i}' for i in self.irange()]
        self._half_space_chart = self.chart(' '.join(xs) + ':(0,oo)')
        self._ball_chart = None
        self._half_space_to_ball = None

        g = self.metric('g')
        eX = self._half_space_chart.frame()
        for i in self.irange():
            g[eX,i,i] = metric_constant^2 / SR.var(xs[-1])^2
        
        self._metric_constant = metric_constant

    def half_space_chart(self):
        return self._half_space_chart        
    
    def ball_chart(self):
        if self._ball_chart == None:
            xs = [f'x{i}' for i in self.irange()]
            ys = [f'y{i}' for i in self.irange()]
            self._ball_chart = self.chart(' '.join(ys), coord_restrictions=(lambda *args: sum(y^2 for y in args) < 1))

            x_denom = sum(SR.var(x)^2 for x in xs[:-1]) + (1 + SR.var(xs[-1]))^2
            x_to_y = [2*SR.var(x) / x_denom for x in xs[:-1]] + [(sum(SR.var(x)^2 for x in xs) - 1) / x_denom]

            y_denom = sum(SR.var(y)^2 for y in ys[:-1]) + (1 - SR.var(ys[-1]))^2
            y_to_x = [2*SR.var(y) / y_denom for y in ys[:-1]] + [(1 + sum(-SR.var(y)^2 for y in ys)) / y_denom]

            self._half_space_to_ball = self._half_space_chart.transition_map(self._ball_chart, x_to_y)
            self._half_space_to_ball.set_inverse(*y_to_x)

        return self._ball_chart

    def half_space_to_ball(self):
        self.ball_chart()
        return self._half_space_to_ball

    def ball_to_half_space(self):
        self.ball_chart()
        return self._half_space_to_ball.inverse()

    def cosmological_constant(self):
        n = self.dim()
        return -(n - 1) * (n - 2) / (2 * self._metric_constant^2)

