#!/usr/bin/env python

import mdtraj as md
import matplotlib.pyplot as plt

t = md.load('npt.trr', top='berendsen.gro')

d = md.density(t)

plt.plot(d)
plt.show()
