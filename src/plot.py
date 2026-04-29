import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    my_data = np.genfromtxt('data/super100.csv', delimiter=',')
    print(my_data[0])
    print(my_data[1])

    print(np.gradient(my_data[1]))

    plt.plot(my_data[0][0:3000], my_data[1][0:3000], label="β = 1/kT",alpha=0.8)
    plt.grid(alpha=0.2)
    plt.title("Vortex Density Phase Transition")

    plt.axvline(x=1/0.637, color="g", linestyle=":", alpha=0.5, label="Ideal theory")
    plt.axvline(x=1/0.893, color="r", linestyle=":", alpha=0.5, label="Raghav G. Jha")


    plt.ylabel("Vortex Density")
    plt.xlabel("Beta")
    plt.xscale("log")
    plt.legend()
    plt.show()

