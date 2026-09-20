"""Build the labelled ROI volume + colour table that BrainNet Viewer renders.

Takes the AAL3 atlas and collapses it into a handful of schematic regions,
renumbered 1..N. BrainNet Viewer draws one isosurface per label value and
colours row i of the colour table to the i-th drawn label, so the label
numbering here and the row order of the colour table have to agree -- both are
driven by the order of the REGIONS lists below.

Two variants are written, because "striatum" and "caudate/putamen/NAcc" are the
same voxels at different granularity and cannot both be drawn:

    split     -- caudate, putamen and nucleus accumbens as three regions
    striatum  -- the three merged into one

Run:  python3 make_rois.py
"""

import numpy as np
import nibabel as nib
from scipy.ndimage import binary_erosion
from pathlib import Path

HERE = Path(__file__).resolve().parent
ATLAS_DIR = HERE / "atlas" / "AAL3"
OUT_DIR = HERE / "rois"

# '1mm' gives visibly smoother blobs; '2mm' renders several times faster in
# MATLAB. Switch here and re-run -- nothing else needs to change.
RESOLUTION = "1mm"

# AAL3 is a partition, so neighbouring regions share a face and render as one
# continuous coloured mass. Shrinking each region by a shell opens a visible
# gap between them. 0 disables. Purely cosmetic -- it makes every region
# slightly smaller than its true anatomical extent, which is fine for a
# schematic but would not be for a quantitative figure.
SEPARATION_MM = 1.0

ATLAS = {"1mm": "AAL3v1_1mm.nii", "2mm": "AAL3v1.nii"}[RESOLUTION]

# (display name, AAL3 label indices to merge, hex colour)
#
# dlPFC, mPFC and vmPFC are functional names with no AAL3 label of their own;
# the anatomical proxies used here are noted alongside. Name the atlas and
# these proxies in the figure caption.
REGIONS_SPLIT = [
    ("dlPFC",   [5, 6],       "#D62728"),   # Frontal_Mid_2 (middle frontal gyrus)
    ("mPFC",    [19, 20],     "#2CA02C"),   # Frontal_Sup_Medial
    ("vmPFC",   [21, 22],     "#8CC63F"),   # Frontal_Med_Orb
    ("dACC",    [155, 156],   "#1F77B4"),   # ACC_sup (supracallosal)
    ("pACC",    [153, 154],   "#5DADE2"),   # ACC_pre (pregenual)
    ("sgACC",   [151, 152],   "#9467BD"),   # ACC_sub (subgenual)
    ("Caudate", [75, 76],     "#FF9900"),
    ("Putamen", [77, 78],     "#E2571E"),
    ("NAcc",    [157, 158],   "#FFD24D"),
]

REGIONS_STRIATUM = REGIONS_SPLIT[:6] + [
    ("Striatum", [75, 76, 77, 78, 157, 158], "#F0862E"),
]

# Striatum alone, no cortex. Colours are carried over from REGIONS_SPLIT so the
# same structure reads the same colour across every figure panel.
REGIONS_STRIATAL = REGIONS_SPLIT[6:]

VARIANTS = {
    "split": REGIONS_SPLIT,
    "striatum": REGIONS_STRIATUM,
    "striatal": REGIONS_STRIATAL,
}


def hex_to_rgb01(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def shrink(mask, n_iter, name):
    """Erode by n_iter voxels, backing off if the region would be gutted."""
    if n_iter < 1:
        return mask
    before = int(mask.sum())
    for k in range(n_iter, 0, -1):
        eroded = binary_erosion(mask, iterations=k)
        kept = int(eroded.sum())
        # Small structures (NAcc, sgACC) can lose most of their volume to a
        # gap that barely dents the big cortical slabs.
        if kept > 0.4 * before:
            if k < n_iter:
                print(f"      ({name}: eroded {k} not {n_iter} voxels to keep volume)")
            return eroded
    print(f"      ({name}: too small to erode, left intact)")
    return mask


def build(variant, regions, atlas_img, atlas_data, erode_iter):
    out = np.zeros(atlas_data.shape, dtype=np.uint8)

    for new_idx, (name, src_labels, _) in enumerate(regions, start=1):
        mask = np.isin(atlas_data, src_labels)
        n_vox = int(mask.sum())
        if n_vox == 0:
            raise ValueError(
                f"{name}: AAL3 labels {src_labels} matched no voxels -- "
                f"check the indices against {ATLAS_DIR / (ATLAS + '.txt')}"
            )
        mask = shrink(mask, erode_iter, name)
        # Later regions must not silently eat earlier ones. AAL3 is a
        # partition so this should never fire, but merging groups by hand is
        # exactly where an overlap would sneak in.
        clash = np.logical_and(mask, out > 0)
        if clash.any():
            raise ValueError(f"{name} overlaps an earlier region in {clash.sum()} voxels")
        out[mask] = new_idx
        print(f"  {new_idx:2d}  {name:<9s} {int(mask.sum()):>7d} voxels "
              f"(of {n_vox})  {src_labels}")

    stem = OUT_DIR / f"schematic_rois_{variant}"

    img = nib.Nifti1Image(out, atlas_img.affine)
    img.set_data_dtype(np.uint8)
    nib.save(img, str(stem) + ".nii")

    colours = np.array([hex_to_rgb01(c) for _, _, c in regions])
    np.savetxt(str(stem) + "_colors.txt", colours, fmt="%.6f", delimiter="  ")

    with open(str(stem) + "_labels.txt", "w") as f:
        f.write("# index\tname\thex\tAAL3_labels\n")
        for i, (name, src, hexc) in enumerate(regions, start=1):
            f.write(f"{i}\t{name}\t{hexc}\t{','.join(map(str, src))}\n")

    return stem, colours


def legend(regions, colours, path):
    """Swatch strip to paste next to the render when assembling the figure."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available -- skipping legend)")
        return

    fig, ax = plt.subplots(figsize=(2.2, 0.32 * len(regions)))
    for i, ((name, _, _), rgb) in enumerate(zip(regions, colours)):
        y = len(regions) - i - 1
        ax.add_patch(plt.Rectangle((0, y + 0.15), 0.5, 0.7, color=rgb))
        ax.text(0.68, y + 0.5, name, va="center", fontsize=10)
    ax.set_xlim(0, 2.2)
    ax.set_ylim(0, len(regions))
    ax.axis("off")
    fig.savefig(path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  legend -> {path.name}")


def main():
    OUT_DIR.mkdir(exist_ok=True)

    atlas_path = ATLAS_DIR / ATLAS
    if not atlas_path.exists():
        raise SystemExit(f"AAL3 atlas not found at {atlas_path}")

    atlas_img = nib.load(str(atlas_path))
    atlas_data = np.asarray(atlas_img.dataobj).astype(np.int32)

    vox_mm = float(np.abs(atlas_img.header.get_zooms()[0]))
    erode_iter = int(round(SEPARATION_MM / vox_mm))
    print(f"atlas: {ATLAS}  shape={atlas_data.shape}  voxel={vox_mm:g}mm")
    print(f"separation: {SEPARATION_MM:g}mm -> erode {erode_iter} voxel(s)")

    for variant, regions in VARIANTS.items():
        print(f"\n[{variant}]")
        stem, colours = build(variant, regions, atlas_img, atlas_data, erode_iter)
        legend(regions, colours, Path(str(stem) + "_legend.png"))
        print(f"  -> {stem.name}.nii  (+ _colors.txt, _labels.txt)")


if __name__ == "__main__":
    main()
