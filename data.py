
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

 # ------------------- Read parameters from input.in -------------------
def read_params(filename):
    params = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Remove inline comments
            line = line.split("#")[0].strip()
            parts = line.split()
            if len(parts) >= 2:
                key = parts[0]
                val = parts[1]
                # Try to parse numeric values, otherwise keep as string
                try:
                    params[key] = float(val)
                except ValueError:
                    params[key] = val
    return params


if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = "input-16.in"
params = read_params(input_file)

# T_target may be numeric in the input file
T_target = float(params.get("T_target", 300.0))  # Default 300K

# Determine output file name from input params (fall back to output.log)
output_param = params.get("output", None)
if output_param is None:
    output_param = "output.log"

# If output_param is not an absolute path, resolve relative to the input file directory
input_dir = os.path.dirname(input_file) if os.path.dirname(input_file) != "" else "."
output_path = output_param if os.path.isabs(output_param) else os.path.join(input_dir, output_param)

print(f"Using input file: {input_file}")
print(f"Using output file: {output_path}")

 # ------------------- Read simulation results -------------------
try:
    data = np.loadtxt(output_path, comments="#", skiprows=1)
except Exception as e:
    print(f"Error reading output file '{output_path}': {e}")
    sys.exit(1)

steps = data[:, 0]
epot = data[:, 1]
ekin = data[:, 2]
etot = data[:, 3]
temp = data[:, 4]
avg_temp = np.mean(temp)
rms_temp = np.sqrt(np.mean((temp - T_target)**2))
avg_ekin = np.mean(ekin)
rms_ekin = np.sqrt(np.mean((ekin - np.mean(ekin))**2))

print(f"Average temperature: {avg_temp:.5f}")
print(f"Temperature RMS deviation from target: {rms_temp:.5f}")
print(f"Average kinetic energy: {avg_ekin:.5f}")
print(f"Kinetic energy RMS deviation: {rms_ekin:.5f}")


 # ------------------- Plotting -------------------
fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
 # Kinetic energy curve
axes[0].plot(steps, ekin, label="E_kin")
axes[0].set_ylabel("Energy")
axes[0].set_title("Energy vs Steps")
axes[0].legend(loc='upper right')  # Legend in the upper right corner
axes[0].grid(True)

 # Annotate average and RMS (top left)
axes[0].text(0.02, 0.95, f"Avg E_kin = {avg_ekin:.5f}\nRMS = {rms_ekin:.5f}",
             transform=axes[0].transAxes,
             verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6))

 # Temperature curve
axes[1].plot(steps, temp, label="Temperature", color="tab:red")
axes[1].axhline(T_target, color="tab:blue", linestyle="--", label=f"Target Temp = {T_target}")
axes[1].set_xlabel("Step")
axes[1].set_ylabel("Temperature")
axes[1].set_title("Temperature vs Steps")
axes[1].legend(loc='upper right')  # Legend in the upper right corner
axes[1].grid(True)

 # Annotate average temperature and RMS (top left)

axes[1].text(0.02, 0.95, f"Avg Temp = {avg_temp:.5f} \nRMSE = {rms_temp:.5f} ",
             transform=axes[1].transAxes,
             verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6))

plt.tight_layout()
plt.savefig("md_energy_temp.png", dpi=300)
