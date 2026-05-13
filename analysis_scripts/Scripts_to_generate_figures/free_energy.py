from PIL import Image, ImageOps
import matplotlib.pyplot as plt

# ---- INPUT FILES ----
files = [
    r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL/Fig3A.png",
    r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL/Fig3B.png",
    r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL/Fig3C.png",
    r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL/Fig3D.png"
]

# ---- PARAMETERS ----
scale_factor = 1.8   # increases effective font size
target_size = (1200, 1200)  # force consistent square panels

processed_images = []

for f in files:
    img = Image.open(f).convert("RGB")

    # Trim whitespace
    img = ImageOps.crop(img, border=10)

    # Resize (increase readability)
    new_size = (int(img.width * scale_factor), int(img.height * scale_factor))
    img = img.resize(new_size, Image.LANCZOS)

    # Force same size for all panels
    img = ImageOps.fit(img, target_size, Image.LANCZOS)

    processed_images.append(img)

# ---- CREATE FIGURE ----
fig, axes = plt.subplots(2, 2, figsize=(7, 7))  # ~double column JCIM

labels = ['A', 'B', 'C', 'D']

for ax, img, label in zip(axes.flatten(), processed_images, labels):
    ax.imshow(img)
    ax.axis('off')

    # Panel labels (JCIM style)
    ax.text(
        0.02, 0.95, label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight='bold',
        va='top',
        ha='left'
    )

# Tight layout (removes gaps)
plt.subplots_adjust(wspace=0.02, hspace=0.02)

# ---- SAVE ----
output_path = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL/Fig3_JCIM.png"
plt.savefig(output_path, dpi=600, bbox_inches='tight')

print("Final JCIM-ready figure saved at:")
print(output_path)

