import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

def crop_whitespace(img, threshold=245):
    if img.ndim == 3:
        gray = img.mean(axis=2)
    else:
        gray = img
    mask = gray < threshold
    coords = np.argwhere(mask)
    if coords.size == 0:
        return img
    y0, x0 = coords.min(axis=0)
    y1, x1 = coords.max(axis=0)
    return img[y0:y1+1, x0:x1+1]


# Paths
imgA_path = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Introduction\Presentation1\Slide1.PNG"
imgB_path = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Introduction\Presentation1\Slide2.PNG"

# Load & crop
imgA = crop_whitespace(mpimg.imread(imgA_path))
imgB = crop_whitespace(mpimg.imread(imgB_path))

hA, wA = imgA.shape[:2]
hB, wB = imgB.shape[:2]

fig_width = 20  # inches
scale     = fig_width / max(wA, wB)

panelA_h = hA * scale
panelB_h = hB * scale
total_h   = panelA_h + panelB_h

fig = plt.figure(figsize=(fig_width, total_h))

fracA = panelA_h / total_h
fracB = panelB_h / total_h

# Pixel offset from top-left corner of each image
PAD_X =  30   # pixels from left edge  ← tune this
PAD_Y =  30   # pixels from top edge   ← tune this

# Panel A — top
ax1 = fig.add_axes([0, fracB, 1, fracA])
ax1.imshow(imgA)
ax1.axis("off")
ax1.text(
    PAD_X, PAD_Y, "(a)",
    fontfamily="Arial",
    fontsize=15,
    color="black",
    ha="left", va="top",
)

# Panel B — bottom
ax2 = fig.add_axes([0, 0, 1, fracB])
ax2.imshow(imgB)
ax2.axis("off")
ax2.text(
    PAD_X, PAD_Y, "(b)",
    fontfamily="Arial",
    fontsize=15,
    color="black",
    ha="left", va="top",
)

plt.savefig(
    r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Introduction\Merged_AB.png",
    dpi=600,
    bbox_inches="tight",
    pad_inches=0,
)

plt.show()