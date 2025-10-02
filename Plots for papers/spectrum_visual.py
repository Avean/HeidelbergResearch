import numpy as np
import matplotlib.pyplot as plt

# Settings
N = 10  # 10 unique + 10 conjugates = 20 total

# Real parts: centered near -0.3, allow some positives
real_parts = np.random.normal(loc=-0.75, scale=1, size=N)
# real_parts = np.clip(real_parts, -1.5, 1.5)

# Imaginary parts: small values to align near horizontal line
imag_parts = np.random.normal(loc=0.0, scale=0.6, size=N)

# Set a few imaginary parts to exactly zero (real spectrum points)
zero_imag_indices = np.random.choice(N, size=2, replace=False)
imag_parts[zero_imag_indices] = 0.0

# Create unique complex points
spectrum_half = real_parts + 1j * imag_parts

# Add complex conjugates for symmetry
spectrum = np.concatenate([spectrum_half, np.conj(spectrum_half)])

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(0, color='gray', linewidth=1)
ax.axvline(0, color='gray', linewidth=1)

ax.scatter(spectrum.real, spectrum.imag, color='blue', s=80)
ax.plot([-0.1, 0.5], [0, 0], color='blue', linewidth=5, linestyle='-', label='Ref Line')

# Formatting
ax.set_aspect('equal')
ax.grid(True, linestyle='--', linewidth=0.5)
ax.tick_params(labelsize=14)
ax.set_xlabel('Real Part', fontsize=16)
ax.set_ylabel('Imaginary Part', fontsize=16)
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.6, 1.6)
plt.tight_layout()

# Save plot
output_path = "spectrum_horizontal_real_bias.png"
plt.savefig(output_path, dpi=300)
print(f"Plot saved as: {output_path}")
