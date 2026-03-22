import numpy as np
import matplotlib.pyplot as plt


def load_data(fname):
    data = np.loadtxt(fname)
    idx  = np.argsort(data[:, 0])
    data = data[idx]
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]


def find_Tc(temps, Cv, chi):
    Tc_Cv  = temps[np.argmax(Cv)]
    Tc_chi = temps[np.argmax(chi)]
    print(f"Tc from Cv  : {Tc_Cv:.4f}")
    print(f"Tc from chi : {Tc_chi:.4f}")
    print(f"Tc theory   : 2.2691")
    return Tc_Cv, Tc_chi


def plot_quantities(temps, E, mag, Cv, chi):
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    axs[0,0].plot(temps, E,   'o-', color='blue')
    axs[0,0].set(xlabel='T', ylabel='E/N²',   title='Energy vs T')

    axs[0,1].plot(temps, mag, 'o-', color='red')
    axs[0,1].set(xlabel='T', ylabel='|M|/N²', title='Magnetization vs T')

    axs[1,0].plot(temps, Cv,  'o-', color='green')
    axs[1,0].set(xlabel='T', ylabel='Cv',     title='Specific Heat vs T')

    axs[1,1].plot(temps, chi, 'o-', color='purple')
    axs[1,1].set(xlabel='T', ylabel='χ',      title='Susceptibility vs T')

    for ax in axs.flat:
        ax.grid(True)

    plt.tight_layout()
    plt.savefig('ising_thermo.png', dpi=300)
    plt.show()


temps, E, mag, Cv, chi = load_data("ising_results_ferro_PB.txt")
Tc_Cv, Tc_chi = find_Tc(temps, Cv, chi)
plot_quantities(temps, E, mag, Cv, chi)