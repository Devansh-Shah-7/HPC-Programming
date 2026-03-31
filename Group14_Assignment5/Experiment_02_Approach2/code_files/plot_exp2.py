import csv
import matplotlib.pyplot as plt
import numpy as np

# Read data
data = []
with open('exp2_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append({
            'Approach': int(row['Approach']),
            'Threads': int(row['Threads']),
            'Total_Time': float(row['Total_Time'])
        })

# Separate by approach
def_data = [d for d in data if d['Approach'] == 1]
imm_data = [d for d in data if d['Approach'] == 2]

def_data.sort(key=lambda x: x['Threads'])
imm_data.sort(key=lambda x: x['Threads'])

# Get baseline (1 thread time)
def_baseline = next(d for d in def_data if d['Threads'] == 1)['Total_Time']
imm_baseline = next(d for d in imm_data if d['Threads'] == 1)['Total_Time']

# Calculate speedups
def_threads = [d['Threads'] for d in def_data]
def_speedup = [def_baseline / d['Total_Time'] for d in def_data]
imm_threads = [d['Threads'] for d in imm_data]
imm_speedup = [imm_baseline / d['Total_Time'] for d in imm_data]

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(def_threads, def_speedup, 'o-', linewidth=2, markersize=10, 
        label='Deferred Insertion', color='#2E86AB')
ax.plot(imm_threads, imm_speedup, 's-', linewidth=2, markersize=10, 
        label='Immediate Replacement', color='#A23B72')

# Ideal speedup line
max_threads = max(def_threads + imm_threads)
ax.plot([1, max_threads], [1, max_threads], '--', color='gray', 
        alpha=0.7, label='Ideal Speedup')

ax.set_xlabel('Number of Threads', fontsize=12)
ax.set_ylabel('Speedup', fontsize=12)
ax.set_title('Experiment 2: OpenMP Speedup (14M particles, 500×200 grid)', 
             fontsize=14, fontweight='bold')
ax.set_xticks([1, 2, 4, 8, 16])
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('exp2_speedup.png', dpi=300, bbox_inches='tight')
print("Saved: exp2_speedup.png")
plt.show()