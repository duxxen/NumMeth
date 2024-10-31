from typing import List
import matplotlib.pyplot as plt
import pandas as pd

pathq = '../data/quad/'
paths = '../data/simp/'
bform = pd.read_csv(pathq + 'init.csv', header=0, delimiter=';', index_col=False)

def plotnow(form: pd.DataFrame, keys: List[str], title='', color = 'red', pos = 'vert', single = False):
    if not single:
        n = len(keys)
        match pos:
            case 'vert':
                axscount = (n, 1)
            case 'horiz':
                axscount = (1, n)
            case 'quad':
                axscount = (n, n)
        
        x, y = axscount
        fig, axs = plt.subplots(x, y, figsize = (n * 3, n * 2))
        axs[0].set_title(title)
        for i in range(0, n):
            key = keys[i]
            axs[i].plot(form['x'], form[key], color=color, label=key)
            axs[i].legend()
            axs[i].grid()
        
    else:
        fig, axs = plt.subplots(figsize = (6, 4))
        axs.set_title(title)
        for key in keys:
            axs.plot(form['x'], form[key], label=key)
        axs.legend()
        axs.grid()

if __name__ == '__main__':

    qtcform = pd.read_csv(pathq + 'test_const_4.csv', header=0, delimiter=';', index_col=False)
    qtfform = pd.read_csv(pathq + 'test_func_4.csv', header=0, delimiter=';', index_col=False)
    stcform = pd.read_csv(paths + 'test_const_4.csv', header=0, delimiter=';', index_col=False)
    stfform = pd.read_csv(paths + 'test_func_4.csv', header=0, delimiter=';', index_col=False)
    bform['quhc'] = qtcform['uhc']
    bform['quhf'] = qtfform['uhf']
    bform['suhc'] = stcform['uhc']
    bform['suhf'] = stfform['uhf']
    
    plotnow(bform, ['u(x)', 'du(x)', 'ddu(x)'], title='Искомая функция u(x)')
    plotnow(bform, ['p(x)', 'q(x)'], title='Параметры p(x), q(x)')
    plotnow(bform, ['f(x)', 'fc(x)'], title='Функция правой части f(x)')
    plotnow(bform, ['u(x)', 'quhc', 'quhf'], single=True, title='Расчет с quad()')
    plotnow(bform, ['u(x)', 'suhc', 'suhf'], single=True, title='Расчет с simp()')
    plt.show()
    

