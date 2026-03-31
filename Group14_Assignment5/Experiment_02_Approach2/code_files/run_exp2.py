import subprocess
import re

# Header
with open('exp2_results.csv', 'w') as f:
    f.write("NX,NY,NUM_Points,Approach,Threads,Interp_Time,Mover_Time,Total_Time\n")

# Read main2.cpp
with open('main2.cpp', 'r') as f:
    original_code = f.read()

for approach in [1, 2]:
    for threads in [1, 2, 4, 8, 16]:
        # Modify code
        code = original_code
        code = re.sub(r'int approach = \d+', f'int approach = {approach}', code)
        code = re.sub(r'int num_threads = \d+', f'int num_threads = {threads}', code)
        
        with open('main2.cpp', 'w') as f:
            f.write(code)
        
        # Compile
        result = subprocess.run(['g++', 'main2.cpp', 'init.cpp', 'utils.cpp', 
                                '-lm', '-fopenmp', '-O2', '-o', 'main2'], 
                               capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Compile error: {result.stderr}")
            continue
        
        # Run
        result = subprocess.run(['./main2'], capture_output=True, text=True)
        output = result.stdout
        
        # Extract CSV
        for line in output.split('\n'):
            if line.startswith('CSV:'):
                csv_line = line.replace('CSV: ', '')
                with open('exp2_results.csv', 'a') as f:
                    f.write(csv_line + '\n')
                print(f"Done: approach={approach}, threads={threads} -> {csv_line}")
                break

print("\nAll done! Check exp2_results.csv")