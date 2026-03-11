import pandas as pd
import matplotlib.pyplot as plt

data=pd.read_csv("exp3_serial_labpc.csv")

plt.plot(data["Iter"],data["Interp"],label="Interp")
plt.plot(data["Iter"],data["Mover"],label="Mover")
plt.plot(data["Iter"],data["Total"],label="Total")

plt.xlabel("Iteration")
plt.ylabel("Time")
plt.legend()
plt.title("Exp-03 Serial")

plt.savefig("exp3_serial_plot.png")
plt.show()
