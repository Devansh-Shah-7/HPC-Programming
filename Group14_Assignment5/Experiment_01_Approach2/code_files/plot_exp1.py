import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data
df = pd.read_csv('exp1_results.csv')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 11

# Create 3 separate plots for each grid size
grids = [(250, 100), (500, 200), (1000, 400)]

for nx, ny in grids:
    fig, ax = plt.subplots()
    
    # Filter data for this grid
    grid_data = df[(df['NX'] == nx) & (df['NY'] == ny)]
    
    # Separate by approach
    deferred = grid_data[grid_data['Approach'] == 1]
    immediate = grid_data[grid_data['Approach'] == 2]
    
    # Plot both approaches
    ax.loglog(deferred['NUM_Points'], deferred['Total_Time'], 
              'o-', linewidth=2, markersize=8, label='Deferred Insertion', color='#2E86AB')
    ax.loglog(immediate['NUM_Points'], immediate['Total_Time'], 
              's-', linewidth=2, markersize=8, label='Immediate Replacement', color='#A23B72')
    
    # Add reference line for O(N) scaling
    x_ref = np.array([1e2, 1e8])
    # Scale to match the middle data point
    mid_point = deferred[deferred['NUM_Points'] == 1e6]['Total_Time'].values[0]
    y_ref = mid_point * (x_ref / 1e6)  # Linear scaling from 1e6 point
    ax.loglog(x_ref, y_ref, '--', color='gray', alpha=0.7, label='O(N) reference')
    
    # Formatting
    ax.set_xlabel('Number of Particles (N)', fontsize=12)
    ax.set_ylabel('Total Time (seconds)', fontsize=12)
    ax.set_title(f'Experiment 1: Serial Scaling ({nx}×{ny} Grid)', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', frameon=True)
    ax.grid(True, which="both", ls="-", alpha=0.3)
    
    # Set axis limits with some padding
    ax.set_xlim(50, 2e8)
    
    plt.tight_layout()
    plt.savefig(f'exp1_grid_{nx}x{ny}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: exp1_grid_{nx}x{ny}.png")
    plt.close()
    plt.show()
print("\nAll 3 plots generated successfully!")