import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.gridspec import GridSpec
import os

# ---- Configuration ----
DPI_SAVE = 300  # Final output DPI

# ---- Image paths ----
image_paths = [
    r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER\1_LFER_distance_colored_human.png",   # A
    r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER\1_LFER_distance_colored_mouse.png",   # B
    r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER\3_Distance_vs_effect_human.png",      # C
    r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER\3_Distance_vs_effect.png",            # D
    r"D:\PhD_Thesis\GPX6\analysis\comparison_dg_star_mouse_vs_human.png",            # E (contains BOTH mouse & human)
    r"D:\PhD_Thesis\GPX6\analysis\comparison_dg0_mouse_vs_human.png",                # F (contains BOTH mouse & human)
]

labels = ["A", "B", "C", "D", "E", "F"]

print("="*80)
print("COMBINED FIGURE - PROPER HEATMAP HANDLING")
print("="*80)

# ---- Load and process images ----
images = []
for i, p in enumerate(image_paths):
    if not os.path.exists(p):
        print(f"WARNING: File not found: {p}")
        continue
        
    file_size_mb = os.path.getsize(p) / (1024 * 1024)
    print(f"\n{labels[i]}: {os.path.basename(p)}")
    print(f"  File size: {file_size_mb:.1f} MB")
    
    img = Image.open(p)
    w, h = img.size
    print(f"  Dimensions: {w} x {h} pixels")
    
    # For heatmaps (E, F), if they're huge, downsample intelligently
    if i >= 4 and (w > 8000 or h > 8000):
        max_dim = 5000  # Slightly larger max for heatmaps
        if w > max_dim or h > max_dim:
            scale = min(max_dim / w, max_dim / h)
            new_w = int(w * scale)
            new_h = int(h * scale)
            print(f"  Downsampling from {w}x{h} to {new_w}x{new_h}")
            img = img.resize((new_w, new_h), Image.Resampling.LANCZOS)
    
    # Keep general plots at original size or slight reduction
    elif w > 6000 or h > 6000:
        scale = min(6000 / w, 6000 / h)
        new_w = int(w * scale)
        new_h = int(h * scale)
        print(f"  Resizing from {w}x{h} to {new_w}x{new_h}")
        img = img.resize((new_w, new_h), Image.Resampling.LANCZOS)
    
    images.append(img)
    print(f"  Loaded: {img.size[0]} x {img.size[1]}")

if len(images) != 6:
    print(f"\nERROR: Only loaded {len(images)} of 6 images!")
    exit(1)

print("\n" + "="*80)
print("Creating combined figure...")
print("="*80)
print("\nNOTE: E and F are FULL WIDTH heatmaps (each contains mouse+human)")

# ---- Create figure with proper layout ----
# Rows 1-2: Two columns each (A-B, C-D)
# Row 3: Single column spanning full width for each heatmap (E, then F below)

fig = plt.figure(figsize=(24, 28), dpi=100)  # Taller to accommodate stacked heatmaps

# Use GridSpec with 4 rows instead of 3
gs = GridSpec(
    4, 2,  # 4 rows, 2 columns
    height_ratios=[1.0, 1.0, 1.4, 1.4],  # Rows 3 and 4 for heatmaps
    hspace=0.04,
    wspace=0.03
)

# Top row: A, B
ax_A = fig.add_subplot(gs[0, 0])
ax_B = fig.add_subplot(gs[0, 1])

# Second row: C, D
ax_C = fig.add_subplot(gs[1, 0])
ax_D = fig.add_subplot(gs[1, 1])

# Third row: E spans full width
ax_E = fig.add_subplot(gs[2, :])  # Spans both columns

# Fourth row: F spans full width
ax_F = fig.add_subplot(gs[3, :])  # Spans both columns

# Arrange axes in order
axes = [ax_A, ax_B, ax_C, ax_D, ax_E, ax_F]

# Plot images
for ax, img, lab in zip(axes, images, labels):
    ax.imshow(img, aspect='auto')
    ax.axis("off")
    
    # Add panel labels
    ax.text(
        0.015, 0.98, lab,
        transform=ax.transAxes,
        fontsize=24,
        fontweight="bold",
        va="top",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="black", linewidth=2.5, pad=6)
    )

# Adjust layout
plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

# ---- Save ----
out = r"D:\PhD_Thesis\GPX6\analysis\GPX6_LFER_Combined_FINAL.png"
print(f"\nSaving to: {out}")
print(f"Output DPI: {DPI_SAVE}")

plt.savefig(out, dpi=DPI_SAVE, bbox_inches="tight", facecolor='white', pad_inches=0.05)
plt.close()

output_size_mb = os.path.getsize(out) / (1024 * 1024)
print(f"\n{'='*80}")
print(f"SUCCESS!")
print(f"{'='*80}")
print(f"Output file: {out}")
print(f"File size: {output_size_mb:.1f} MB")
print(f"\nLayout:")
print(f"  Row 1: A | B")
print(f"  Row 2: C | D")
print(f"  Row 3: E (full width - contains mouse & human ΔG*)")
print(f"  Row 4: F (full width - contains mouse & human ΔG0)")
print(f"{'='*80}")