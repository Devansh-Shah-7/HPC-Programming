import pandas as pd

s = pd.read_csv("exp3_serial_labpc.csv")
p = pd.read_csv("exp3_parallel_labpc.csv")

speedup = s["Mover"] / p["Mover"]

print("\nSpeedup per iteration:\n")
print(speedup)

print("\nAverage Speedup = ", speedup.mean())
