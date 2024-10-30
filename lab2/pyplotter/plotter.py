from typing import List
import matplotlib.pyplot as plt
import pandas as pd

path = '../data/'
bform = pd.read_csv(path + 'init.csv', header=0, delimiter=';', index_col=False)

def plotnow(form: pd.DataFrame, keys: List[str], color = 'red', pos = 'vert', single = False):
    n = len(keys)
    match pos:
        case 'vert':
            axscount = (n, 1)
        case 'horiz':
            axscount = (1, n)
        case 'quad':
            axscount = (n, n)
    
    x, y = axscount
    fig, axs = plt.subplots(x, y, figsize = (10, 8))
    for i in range(0, n):
        key = keys[i]
        axs[i].plot(form['x'], form[key], color=color, label=key)
        axs[i].legend()
        axs[i].grid()
    
    plt.show()

if __name__ == '__main__':
    plotnow(bform, ['u(x)', 'du(x)', 'ddu(x)'])
    plotnow(bform, ['p(x)', 'q(x)'])
    plotnow(bform, ['f(x)', 'fc(x)'])
