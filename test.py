import hsgen
import matplotlib.pyplot as plt

X = hsgen.generate_points_on_circle(1000)

plt.scatter([x[0] for x in X], [x[1] for x in X])
plt.axis('equal')

plt.show()
