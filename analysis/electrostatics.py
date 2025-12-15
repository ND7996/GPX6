from pymol import cmd

# ==============================================================================
# STEP 1: LOAD YOUR STRUCTURE FIRST!
# ==============================================================================
# Uncomment and modify this line with your PDB file:
# cmd.load("your_structure.pdb")

# OR if already loaded, just run this script with: @script.py

# ==============================================================================
# RESIDUE DEFINITIONS
# ==============================================================================
active_site_resi = ["49", "83"]
hotspot_resi = ["47", "52", "87", "99", "143", "144", "177", "181"]

# ==============================================================================
# RESET AND SETUP
# ==============================================================================
cmd.hide("everything", "all")
cmd.show("cartoon", "all")
cmd.color("grey80", "all")
cmd.bg_color("white")

# ==============================================================================
# CREATE SELECTIONS
# ==============================================================================
cmd.select("active_site", "resi " + "+".join(active_site_resi))
cmd.select("hotspots", "resi " + "+".join(hotspot_resi))

# ==============================================================================
# STYLE ACTIVE SITE (YELLOW)
# ==============================================================================
cmd.show("sticks", "active_site")
cmd.color("yellow", "active_site")
cmd.set("stick_radius", 0.3, "active_site")

# ==============================================================================
# STYLE HOTSPOTS (RED)
# ==============================================================================
cmd.show("sticks", "hotspots")
cmd.color("red", "hotspots")
cmd.set("stick_radius", 0.3, "hotspots")

# ==============================================================================
# LABELS (ONE PER RESIDUE)
# ==============================================================================
cmd.set("label_size", 14)
cmd.set("label_color", "black")
cmd.set("label_bg_color", "white")
cmd.label("active_site and name CA", "'%s%s' % (resn, resi)")
cmd.label("hotspots and name CA", "'%s%s' % (resn, resi)")

# ==============================================================================
# DISTANCES
# ==============================================================================
cmd.distance("dist", "hotspots and name CA", "active_site and name CA", cutoff=10)
cmd.hide("labels", "dist")
cmd.set("dash_color", "grey50")
cmd.set("dash_width", 2)

# ==============================================================================
# VIEW - VERY IMPORTANT!
# ==============================================================================
cmd.zoom("all")
cmd.center("all")
cmd.deselect()

print("Done! If you don't see anything, press 'A' to auto-zoom or type: zoom all")