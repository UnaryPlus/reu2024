load('metric_utils.sage')

a = var('a', domain='positive')
t = var('t')

H2 = HyperbolicSpace(2, metric_constant=a)
X = H2.half_space_chart()
g = H2.metric()

# Einstein tensor
G = einstein(g)
assert(G == 0)

# Calculate length of curve
c = H2.curve({X: (1, t)}, (t, 1, 2), name='c')
length_c = arc_length(c, g)
assert(length_c == a * log(2))