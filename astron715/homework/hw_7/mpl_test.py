import numpy as np
import matplotlib.pyplot as plt

plt.close(); plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(range(10), range(10))
testing = True
if testing:
# Mess things up
    fig.canvas.draw()
    labels = [item.get_text().zfill(10) for item in ax.get_xticklabels()]
    ax.set_xticklabels(labels)

fig.canvas.draw()


bbox1 = ax.get_xticklabels()[0].get_window_extent()
bbox2 = ax.get_xticklabels()[1].get_window_extent()

print('overlap before', bbox1.overlaps(bbox2))

bbox2.set_points(np.array([[-np.Inf, -np.Inf], [-np.Inf, -np.Inf]]))

print('overlap after', bbox1.overlaps(bbox2))

bbox2.set_points(np.array([[-np.Inf, -np.Inf], [-np.Inf, -np.Inf]]))

print('overlap after', bbox1.overlaps(bbox2))

import myplotting as myplt
reload(myplt)

ax = myplt.delete_overlapping_xlabels3(fig,ax)

plt.show()


