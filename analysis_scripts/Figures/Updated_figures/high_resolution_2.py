"""
merge_slides_to_svg.py
----------------------
Merges Slide1-Slide4 PNGs into a single SVG (2x2 grid) and saves each
slide individually as its own SVG at 600 DPI.

Requirements:
    pip install Pillow numpy

Usage:
    python merge_slides_to_svg.py
"""

import base64
import io
from pathlib import Path

import numpy as np
from PIL import Image

# ── SETTINGS ──────────────────────────────────────────────────────────────────
INPUT_DIR = Path(
    r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"
    r"\Updated_figures\Updated_figures_V3"
)

SLIDES = [
    INPUT_DIR / "Slide1.PNG",
    INPUT_DIR / "Slide2.PNG",
    INPUT_DIR / "Slide3.PNG",
    INPUT_DIR / "Slide4.PNG",
]

DPI            = 600   # physical DPI declared in SVG dimensions
TRIM_TOLERANCE = 10    # how close to 255 a channel must be to count as white
SAFETY_MARGIN  = 10    # pixels of border kept after trim
GAP_PX         = 60    # gap between panels in merged SVG
# ──────────────────────────────────────────────────────────────────────────────


def flatten_to_rgb(img: Image.Image) -> Image.Image:
    """Flatten any mode to RGB on a white background (handles RGBA transparency)."""
    if img.mode in ("RGBA", "LA") or (img.mode == "P" and "transparency" in img.info):
        bg = Image.new("RGB", img.size, (255, 255, 255))
        rgba = img.convert("RGBA")
        bg.paste(rgba, mask=rgba.split()[3])
        return bg
    return img.convert("RGB")


def trim_whitespace(img: Image.Image, name: str = "") -> Image.Image:
    """Remove near-white borders from all four sides."""
    rgb = flatten_to_rgb(img)
    arr = np.array(rgb, dtype=np.uint16)

    threshold = 255 - TRIM_TOLERANCE
    is_bg = (arr[:, :, 0] >= threshold) & \
            (arr[:, :, 1] >= threshold) & \
            (arr[:, :, 2] >= threshold)

    rows = np.where(~is_bg.all(axis=1))[0]
    cols = np.where(~is_bg.all(axis=0))[0]

    orig_w, orig_h = img.size

    if len(rows) == 0 or len(cols) == 0:
        print(f"  [{name}] WARNING: entirely blank, skipping trim.")
        return img

    top    = max(0,          rows[0]  - SAFETY_MARGIN)
    bottom = min(orig_h - 1, rows[-1] + SAFETY_MARGIN)
    left   = max(0,          cols[0]  - SAFETY_MARGIN)
    right  = min(orig_w - 1, cols[-1] + SAFETY_MARGIN)

    print(f"  [{name}] {orig_w}x{orig_h} -> {right-left+1}x{bottom-top+1} px")
    return flatten_to_rgb(img).crop((left, top, right + 1, bottom + 1))


def img_to_b64(img: Image.Image) -> str:
    """Return base64 PNG data URI."""
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def build_svg(viewbox_w, viewbox_h, mm_w, mm_h, image_tags: list) -> str:
    """
    Build SVG bytes with the XML declaration strictly on the very first line.
    Writes as UTF-8 without BOM to avoid the 'XML declaration not at start' error.
    """
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg"',
        '     xmlns:xlink="http://www.w3.org/1999/xlink"',
        '     width="{:.4f}mm" height="{:.4f}mm"'.format(mm_w, mm_h),
        '     viewBox="0 0 {} {}"'.format(viewbox_w, viewbox_h),
        '     version="1.1">',
    ]
    lines += image_tags
    lines.append('</svg>\n')
    return "\n".join(lines)


def single_svg(img: Image.Image) -> str:
    w, h = img.size
    mm_w = w / DPI * 25.4
    mm_h = h / DPI * 25.4
    uri  = img_to_b64(img)
    tags = [
        '  <image x="0" y="0" width="{}" height="{}"'.format(w, h),
        '         xlink:href="{}"'.format(uri),
        '         image-rendering="optimizeQuality"/>',
    ]
    return build_svg(w, h, mm_w, mm_h, tags)


def merged_2x2_svg(images: list) -> str:
    """Arrange 4 images in a 2x2 grid SVG."""
    assert len(images) == 4, "Need exactly 4 images."

    # Normalise all images to the same width
    max_w = max(img.size[0] for img in images)
    resized = []
    for img in images:
        w, h = img.size
        if w != max_w:
            new_h = int(round(h * max_w / w))
            img = img.resize((max_w, new_h), Image.LANCZOS)
        resized.append(img)

    cell_w  = max_w
    heights = [img.size[1] for img in resized]
    row_h   = [max(heights[0], heights[1]),
               max(heights[2], heights[3])]

    total_w = cell_w * 2 + GAP_PX
    total_h = row_h[0] + row_h[1] + GAP_PX
    mm_w    = total_w / DPI * 25.4
    mm_h    = total_h / DPI * 25.4

    origins = [
        (0,               0),
        (cell_w + GAP_PX, 0),
        (0,               row_h[0] + GAP_PX),
        (cell_w + GAP_PX, row_h[0] + GAP_PX),
    ]

    tags = []
    for idx, (img, (ox, oy)) in enumerate(zip(resized, origins)):
        w, h = img.size
        uri  = img_to_b64(img)
        tags += [
            '  <!-- Slide{} -->'.format(idx + 1),
            '  <image x="{}" y="{}" width="{}" height="{}"'.format(ox, oy, w, h),
            '         xlink:href="{}"'.format(uri),
            '         image-rendering="optimizeQuality"/>',
        ]

    return build_svg(total_w, total_h, mm_w, mm_h, tags)


def write_svg(path: Path, content: str):
    """Write SVG as plain UTF-8 bytes (no BOM) so XML declaration is first."""
    path.write_bytes(content.encode("utf-8"))


def main():
    print("=== Loading & trimming slides ===")
    trimmed = []
    for path in SLIDES:
        if not path.exists():
            raise FileNotFoundError("File not found: {}".format(path))
        img = Image.open(path)
        img = trim_whitespace(img, path.stem)
        trimmed.append(img)

        out = path.with_suffix(".svg")
        write_svg(out, single_svg(img))
        print("  Saved -> {}".format(out.name))

    print("\n=== Building merged 2x2 SVG ===")
    merged = merged_2x2_svg(trimmed)
    out_merged = INPUT_DIR / "merged_2x2.svg"
    write_svg(out_merged, merged)
    print("  Saved -> {}".format(out_merged.name))

    print("\nDone!")


if __name__ == "__main__":
    main()