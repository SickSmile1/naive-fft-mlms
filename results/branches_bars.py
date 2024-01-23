###
# adapted from https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.bar.html
###
import numpy as np
import matplotlib.pyplot as plt


species = ('Naive', 'FFT', 'MLMS')
sex_counts = {
    'Branch Misses %': np.array([1,2.6,7.4]),
    'Branches 10E+9': np.array([29.4,8.8,8.28]),
    'cycles front 10E+9': np.array([0.2,0.04,66.1]),
    'cycles back 10E+9': np.array([43,14.8,8.28]),
    'instructions 10E+8': np.array([33.2,9.22,6.61])
}
width = 0.6  # the width of the bars: can also be len(x) sequence


fig, ax = plt.subplots()
bottom = np.ones(3)

fig.set_dpi(300)

for sex, sex_count in sex_counts.items():
    p = ax.bar(species, sex_count, width, label=sex, bottom=bottom)
    bottom += sex_count

    ax.bar_label(p, label_type='center')

ax.yaxis.set_visible(False)
#ax.set_title('Branch, cycle and instruction comparison')
ax.legend()

plt.tight_layout()
plt.show()