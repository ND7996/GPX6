import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

def crop_white(img, tol=245):
    """Crop near-white borders from image."""
    if img.ndim == 3:
        gray = img.mean(axis=2)
    else:
        gray = img

    mask = gray < tol
    coords = np.argwhere(mask)

    y0, x0 = coords.min(axis=0)
    y1, x1 = coords.max(axis=0) + 1

    return img[y0:y1, x0:x1]

# ================================
# Image paths
# ================================
paths = {
    "A": r"D:\PhD_Thesis\Article-GPX6-EVB\Figures\mousecys_FreeEnergyProfile_AllReplicas.png",
    "B": r"D:\PhD_Thesis\Article-GPX6-EVB\Figures\mousesec_FreeEnergyProfile_AllReplicas.png",
    "C": r"D:\PhD_Thesis\Article-GPX6-EVB\Figures\humancys_FreeEnergyProfile_AllReplicas.png",
    "D": r"D:\PhD_Thesis\Article-GPX6-EVB\Figures\humansec_FreeEnergyProfile_AllReplicas.png",
}

save_path = r"D:\PhD_Thesis\Article-GPX6-EVB\Figures\FreeEnergyProfiles.png"

# ================================
# Load & crop
# ================================
imgs = {k: crop_white(mpimg.imread(v)) for k, v in paths.items()}

# ================================
# Create figure
# ================================
fig, axes = plt.subplots(
    2, 2,
    figsize=(12, 10)   # square-ish → less white
)

order = [["A", "B"],
         ["C", "D"]]

for i in range(2):
    for j in range(2):
        key = order[i][j]
        ax = axes[i, j]
        ax.imshow(imgs[key])
        ax.axis("off")

        # Panel labels
        ax.text(
            -0.03, 1.02, key,
            transform=ax.transAxes,
            fontsize=18,
            fontweight="bold",
            va="bottom",
            ha="right"
        )

# ================================
# Tight spacing
# ================================
plt.subplots_adjust(
    left=0.02,
    right=0.98,
    top=0.97,
    bottom=0.03,
    wspace=0.04,
    hspace=0.06
)

# ================================
# Save
# ================================
plt.savefig(
    save_path,
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.01
)

plt.show()
plt.close()
