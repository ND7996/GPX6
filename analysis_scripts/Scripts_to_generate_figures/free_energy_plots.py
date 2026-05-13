"""
unify_xaxis_png.py
------------------
Unifies the x-axis range of 4 existing PNG figures to -350..500
while keeping the y-axis and all other content intact.

Requirements: pip install Pillow
"""
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

# â”€â”€ Spine coordinates (auto-detected from images) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
LEFT, RIGHT = 263, 1786
TOP, BOTTOM = 206, 1416
IMG_W, IMG_H = 2000, 1607
PLOT_W = RIGHT - LEFT   # 1523 px
PLOT_H = BOTTOM - TOP   # 1210 px

# â”€â”€ Input figures with their current x-axis data ranges â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
INPUT_DIR  = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL"
OUTPUT_DIR = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL/unified_xaxis"

fig_info = [
    ("Fig3A.png", -300, 325),   # humancys
    ("Fig3B.png", -300, 400),   # humansec
    ("Fig3C.png", -300, 400),   # mousecys
    ("Fig3D.png", -300, 500),   # mousesec
]

# â”€â”€ Unified target x range â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
TARGET_XMIN, TARGET_XMAX = -350, 500
px_per_unit = PLOT_W / (TARGET_XMAX - TARGET_XMIN)   # ~1.792 px/unit

# â”€â”€ Tick marks to draw â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
TICKS = [-300, -200, -100, 0, 100, 200, 300, 400, 500]

# â”€â”€ ACS-style font (adjust path for Windows if needed) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# Windows: r"/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONT_PATH  = r"/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
TICK_SIZE  = 44
LABEL_SIZE = 44

try:
    font_tick  = ImageFont.truetype(FONT_PATH, TICK_SIZE)
    font_label = ImageFont.truetype(FONT_PATH, LABEL_SIZE)
except:
    font_tick = font_label = ImageFont.load_default()

os.makedirs(OUTPUT_DIR, exist_ok=True)

for fname, xmin, xmax in fig_info:
    fpath = os.path.join(INPUT_DIR, fname)
    if not os.path.exists(fpath):
        print(f"[SKIP] {fpath}")
        continue

    img = Image.open(fpath).convert('RGB')

    # 1. Extract interior pixels (strictly inside spines)
    interior = np.array(img)[TOP+1:BOTTOM, LEFT+1:RIGHT, :]
    int_h, int_w = interior.shape[:2]

    # 2. Rescale interior to match target px/unit
    new_int_w = int((xmax - xmin) * px_per_unit)
    rescaled  = Image.fromarray(interior).resize((new_int_w, int_h), Image.LANCZOS)

    # 3. Compose into white canvas
    canvas = Image.new('RGB', (int_w, int_h), (255, 255, 255))
    left_offset = int((xmin - TARGET_XMIN) * px_per_unit)
    canvas.paste(rescaled, (left_offset, 0))

    # 4. Copy original (preserves title, y-axis, margins)
    result = img.copy()

    # 5. Paste new interior (does NOT touch spine rows or title)
    result.paste(canvas, (LEFT+1, TOP+1))

    # 6. White-out old x-axis strip
    draw = ImageDraw.Draw(result)
    draw.rectangle([(0, BOTTOM+1), (IMG_W, IMG_H)], fill=(255, 255, 255))

    # 7. New tick marks and labels
    TICK_LEN = 14
    for tick_val in TICKS:
        frac = (tick_val - TARGET_XMIN) / (TARGET_XMAX - TARGET_XMIN)
        tx = int(LEFT + frac * PLOT_W)
        draw.line([(tx, BOTTOM), (tx, BOTTOM + TICK_LEN)], fill=(0,0,0), width=3)
        label_text = str(tick_val)
        bbox = draw.textbbox((0,0), label_text, font=font_tick)
        tw = bbox[2] - bbox[0]
        draw.text((tx - tw//2, BOTTOM + TICK_LEN + 6), label_text,
                  fill=(0,0,0), font=font_tick)

    # 8. X-axis label
    xlabel = "E1-E2 [kcal/mol]"
    bbox   = draw.textbbox((0,0), xlabel, font=font_label)
    xw     = bbox[2] - bbox[0]
    xlabel_y = BOTTOM + TICK_LEN + TICK_SIZE + 28
    draw.text(((LEFT+RIGHT)//2 - xw//2, xlabel_y), xlabel,
              fill=(0,0,0), font=font_label)

    out_path = os.path.join(OUTPUT_DIR, fname)
    result.save(out_path, dpi=(300, 300))
    print(f"[OK] {out_path}")

print("\nDone. All figures saved with unified x-axis âˆ’350 to 500.")

