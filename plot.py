import re
import numpy as np
import matplotlib.pyplot as plt

labels = [ 'Dipole', 'Better DP', 'QD', 'MC' ]

def main():
    f = open('output.dat', 'r')
    lines = [ l for l in f ]
    n = len(lines)
    data = np.ndarray((n, 5))
    for i in range(n):
        data[i,:] = [ float(x) for x in re.split('\s+', lines[i].strip()) ]

    fig, ax = plt.subplots()
    lines = [ None ] * 4
    for i in range(4):
        lines[i], = ax.plot(data[:,0], data[:,i+1], label=labels[i])

    ax.legend(handles=lines, loc=1)
    plt.savefig('figure.png', dpi=150)
    plt.show()

if __name__ == '__main__':
    main()
