import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("exp3_serial_cluster.csv")

plt.plot(data["Iter"], data["Interp"], marker='o', label="Interpolation")
plt.plot(data["Iter"], data["Mover"], marker='o', label="Mover")
plt.plot(data["Iter"], data["Total"], marker='o', label="Total")

plt.xlabel("Iteration")
plt.ylabel("Time (sec)")
plt.title("Experiment-03 Cluster")
plt.legend()

plt.savefig("exp3_cluster_plot.png")
plt.show()
