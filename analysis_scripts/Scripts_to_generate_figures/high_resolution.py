from PIL import Image
import numpy as np
import os

input_dir = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Updated_figures\High_resolution"
output_dir = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Updated_figures\High_resolution"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

TARGET_DPI = 600
FIGURE_WIDTH_INCHES = 7.0
WHITE_THRESHOLD = 240  # pixels brighter than this are treated as "white" (0-255)
PADDING_PERCENT = 0.02  # Add 2% padding to ensure nothing is cut

png_files = [f for f in os.listdir(input_dir) if f.lower().endswith(".png")]

for filename in png_files:
    input_path = os.path.join(input_dir, filename)
    name_no_ext = os.path.splitext(filename)[0]
    
    # Define output paths for both formats
    output_tiff_path = os.path.join(output_dir, f"{name_no_ext}_JCIM.tiff")
    output_png_path = os.path.join(output_dir, f"{name_no_ext}_JCIM.png")

    with Image.open(input_path) as img:
        # Flatten transparency onto white background
        if img.mode in ("RGBA", "P", "LA"):
            bg = Image.new("RGB", img.size, (255, 255, 255))
            if img.mode == "RGBA":
                bg.paste(img, mask=img.split()[3])
            else:
                bg.paste(img.convert("RGB"))
            img = bg
        else:
            img = img.convert("RGB")

        # Convert to numpy array
        arr = np.array(img)
        
        # Find non-white pixels (any channel significantly below white)
        non_white = np.any(arr < WHITE_THRESHOLD, axis=2)
        
        # Manual dilation to include anti-aliased edges
        kernel_size = 3
        non_white_dilated = non_white.copy()
        
        # Simple dilation: expand the mask by kernel_size pixels in all directions
        for i in range(-kernel_size, kernel_size + 1):
            for j in range(-kernel_size, kernel_size + 1):
                if i == 0 and j == 0:
                    continue
                shifted = np.roll(non_white, shift=(i, j), axis=(0, 1))
                non_white_dilated = non_white_dilated | shifted
        
        # Find bounding box of content
        rows = np.any(non_white_dilated, axis=1)
        cols = np.any(non_white_dilated, axis=0)
        
        row_indices = np.where(rows)[0]
        col_indices = np.where(cols)[0]
        
        if len(row_indices) > 0 and len(col_indices) > 0:
            row_min, row_max = row_indices[0], row_indices[-1]
            col_min, col_max = col_indices[0], col_indices[-1]
            
            # Calculate padding (percentage of content size)
            pad_h = max(5, int((row_max - row_min) * PADDING_PERCENT))
            pad_w = max(5, int((col_max - col_min) * PADDING_PERCENT))
            
            # Apply padding
            row_min = max(0, row_min - pad_h)
            row_max = min(arr.shape[0], row_max + pad_h)
            col_min = max(0, col_min - pad_w)
            col_max = min(arr.shape[1], col_max + pad_w)
            
            # Crop with padding
            img = img.crop((col_min, row_min, col_max, row_max))
        
        # Get dimensions after cropping
        orig_w, orig_h = img.size
        
        # Calculate target dimensions (width-based scaling)
        aspect_ratio = orig_h / orig_w
        target_w_px = int(FIGURE_WIDTH_INCHES * TARGET_DPI)
        target_h_px = int(target_w_px * aspect_ratio)
        
        # Resize with high quality
        img_resized = img.resize((target_w_px, target_h_px), Image.Resampling.LANCZOS)
        
        # Ensure RGB mode
        if img_resized.mode != "RGB":
            img_resized = img_resized.convert("RGB")
        
        # Save as TIFF (for publication)
        img_resized.save(
            output_tiff_path,
            format="TIFF",
            dpi=(TARGET_DPI, TARGET_DPI),
            compression="tiff_lzw"
        )
        
        # Save as PNG (for preview/screening)
        img_resized.save(
            output_png_path,
            format="PNG",
            dpi=(TARGET_DPI, TARGET_DPI),
            optimize=True
        )
        
        print(f"âœ“ {filename}")
        print(f"  Cropped: {orig_w}x{orig_h} â†’ Resized: {target_w_px}x{target_h_px}")
        print(f"  Padding: {pad_w}px horizontal, {pad_h}px vertical")
        print(f"  Aspect ratio preserved: {aspect_ratio:.3f}")
        print(f"  Saved TIFF: {os.path.basename(output_tiff_path)}")
        print(f"  Saved PNG:  {os.path.basename(output_png_path)}")
        print()

print(f"\nâœ… Complete! {len(png_files)} images processed")
print(f"Output directory: {output_dir}")
print(f"\nJCIM Publication Requirements:")
print(f"  âœ“ Resolution: {TARGET_DPI} DPI")
print(f"  âœ“ Width: {FIGURE_WIDTH_INCHES}\" ({target_w_px}px at {TARGET_DPI} DPI)")
print(f"  âœ“ Format: TIFF with LZW compression")
print(f"  âœ“ PNG also saved for preview purposes")
print(f"  âœ“ Whitespace removed with {PADDING_PERCENT*100}% padding")
print(f"  âœ“ Content preserved (nothing cut)")

