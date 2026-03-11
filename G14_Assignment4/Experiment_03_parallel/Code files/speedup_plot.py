import pandas as pd
import matplotlib.pyplot as plt

serial = pd.read_csv("exp3_serial_labpc.csv")
parallel = pd.read_csv("exp3_parallel_labpc.csv")

speedup = serial["Mover"] / parallel["Mover"]

plt.plot(serial["Iter"], speedup, marker='o')

plt.xlabel("Iteration")
plt.ylabel("Speedup")
plt.title("Speedup of Parallel Mover")

plt.savefig("speedup_plot.png")
plt.show()
