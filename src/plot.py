import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    my_data = np.genfromtxt('data/Best.csv', delimiter=',')
    print(my_data[0])
    print(my_data[1])

    plt.plot(my_data[0], my_data[1])
    plt.show()