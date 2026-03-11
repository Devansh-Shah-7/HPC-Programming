import pandas as pd
import matplotlib.pyplot as plt

serial = pd.read_csv("exp3_serial_labpc.csv")
parallel = pd.read_csv("exp3_parallel_labpc.csv")

plt.plot(serial["Iter"], serial["Mover"], marker='o', label="Serial Mover")
plt.plot(parallel["Iter"], parallel["Mover"], marker='o', label="Parallel Mover")

plt.xlabel("Iteration")
plt.ylabel("Mover Time (sec)")
plt.title("Serial vs Parallel Mover Time")
plt.legend()

plt.savefig("mover_serial_vs_parallel.png")
plt.show()
