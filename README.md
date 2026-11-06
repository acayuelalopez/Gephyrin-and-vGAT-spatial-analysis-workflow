# Gephyrin and vGAT spatial analysis workflow

This workflow detects **Gephyrin** and **vGAT** puncta in ipsilateral and contralateral cortex annotations using **QuPath**, then quantifies their spatial distribution in concentric rings around the tumor using **Python**. Ring counts are normalized by the analyzed ring area and reported as densities in cells/mm².

## Workflow overview

1. Prepare the tumor, ipsilateral cortex, and contralateral cortex annotations in QuPath.
2. Run the QuPath Groovy scripts for Gephyrin and vGAT on each side.
3. Export the tumor and cortex annotations as GeoJSON.
4. Organize the exported CSV and GeoJSON files into one folder per image.
5. Run `analyze_distance_normalized_ring_fast.py`.
6. Review the per-image ring tables, overlays, plots, PDF reports, and pooled summaries.

## Repository contents

```text
.
├── Gephyrin_analysis_contralateral.groovy
├── Gephyrin_analysis_ipsilateral.groovy
├── vGat_analysis_contralateral.groovy
├── vGat_analysis_ipsilateral.groovy
├── analyze_distance_normalized_ring_fast.py
└── README.md
```

Only the contralateral Groovy scripts may be included in some distributions. The ipsilateral scripts follow the same structure, but select the `Cortex_ipsilateral` class and use `ipsilateral` in the exported filenames.

---

# Part 1: QuPath analysis

## Requirements

- QuPath with a compatible scripting environment.
- QuPath Cellpose extension from BIOP.
- A working Cellpose installation configured in QuPath.
- The Cellpose model used by the scripts:

  ```text
  model_telomer.901241
  ```

- Images with spatial calibration in µm/pixel.
- A QuPath project containing the images to analyze.

## Required annotations

Each image should contain the following classified annotations:

```text
Tumor
Cortex_ipsilateral
Cortex_contralateral
```

The Groovy scripts filter annotations by **PathClass**, not merely by the displayed annotation name. The class spelling must therefore match the value used in the script.

For example, the contralateral scripts select:

```groovy
it.getPathClass()?.toString() == "Cortex_contralateral"
```

The ipsilateral versions should select:

```groovy
it.getPathClass()?.toString() == "Cortex_ipsilateral"
```

## Marker-specific detection settings

| Marker | QuPath channel | Selected annotation | Output label |
|---|---:|---|---|
| Gephyrin | `Channel 2` | `Cortex_contralateral` or `Cortex_ipsilateral` | `Gephyrin` |
| vGAT | `Channel 3` | `Cortex_contralateral` or `Cortex_ipsilateral` | `vGat` |

The supplied scripts use the following Cellpose parameters:

```groovy
diameterNuc = 1
flowThreshold = 0.6
cellprobThreshold = -3.0
measureIntensity()
```

The detection diameter is passed to Cellpose after conversion using the image pixel size:

```groovy
.diameter(diameterNuc / pixelSize)
```

These parameters should be validated visually on representative images before processing the complete dataset.

## QuPath outputs

For each marker and side, the scripts export:

```text
<ImageName>_Gephyrin_contralateral_summary.csv
<ImageName>_Gephyrin_contralateral_single_cell.csv
<ImageName>_vGat_contralateral_summary.csv
<ImageName>_vGat_contralateral_single_cell.csv
```

The ipsilateral versions should produce:

```text
<ImageName>_Gephyrin_ipsilateral_summary.csv
<ImageName>_Gephyrin_ipsilateral_single_cell.csv
<ImageName>_vGat_ipsilateral_summary.csv
<ImageName>_vGat_ipsilateral_single_cell.csv
```

The Python analysis uses the files containing `single_cell` in their names. Marker and side are inferred from the filename. If the side is absent from the filename, the script attempts to infer it from the `Parent` column.

## Running the Groovy scripts

For each project image:

1. Open the image in QuPath.
2. Confirm that the relevant cortex annotations have the correct PathClass.
3. Open **Automate > Show script editor**.
4. Run the Gephyrin script for the required side.
5. Inspect the detections.
6. Run the vGAT script for the same side.
7. Repeat for the other side.
8. Save the QuPath project data.

### Important detection behavior

The vGAT script includes:

```groovy
clearDetections()
```

This removes existing detections before running vGAT. The Gephyrin single-cell CSV must therefore be exported successfully before the vGAT script is executed. If detections need to remain visible in QuPath, remove or comment out this line only after confirming that doing so will not mix old and new detections.

## Exporting annotation geometries

The Python script also needs the geometry of:

- `Tumor`
- `Cortex_ipsilateral`
- `Cortex_contralateral`

Export each annotation from QuPath as **GeoJSON**. The filename must contain the corresponding keyword, for example:

```text
Sample_01_Tumor.geojson
Sample_01_Cortex_ipsilateral.geojson
Sample_01_Cortex_contralateral.geojson
```

The Python script can read `.geojson`, `.json`, `.zip`, and compatible `.txt` geometry files. GeoJSON is the recommended format.

---

# Part 2: Input organization for Python

The Python input root must contain **one subfolder per image**. Do not place all CSV files directly in the root folder.

Recommended structure:

```text
analysis_input/
├── Sample_01/
│   ├── Sample_01_Tumor.geojson
│   ├── Sample_01_Cortex_ipsilateral.geojson
│   ├── Sample_01_Cortex_contralateral.geojson
│   ├── Sample_01_Gephyrin_ipsilateral_single_cell.csv
│   ├── Sample_01_Gephyrin_contralateral_single_cell.csv
│   ├── Sample_01_vGat_ipsilateral_single_cell.csv
│   └── Sample_01_vGat_contralateral_single_cell.csv
├── Sample_02/
│   └── ...
└── Sample_03/
    └── ...
```

Summary CSV files exported by QuPath may remain in these folders, but they are not required for the ring analysis.

## Required columns in single-cell CSV files

The Python script automatically searches for columns containing:

```text
Centroid X
Centroid Y
```

Whenever possible, it prioritizes coordinate columns expressed in µm. The CSV files may be comma-separated or tab-separated.

---

# Part 3: Python analysis

## Python requirements

Recommended Python version:

```text
Python 3.10 or newer
```

Required packages:

```text
numpy
pandas
matplotlib
shapely
scipy
```

`scipy` is used for the optional Mann-Whitney U comparisons between ipsilateral and contralateral measurements. The remaining analysis can run without this statistical output if SciPy is unavailable.

Example installation:

```bash
python -m pip install numpy pandas matplotlib shapely scipy
```

## Basic execution

```bash
python analyze_distance_normalized_ring_fast.py \
  --root-dir "/path/to/analysis_input" \
  --output-dir "/path/to/analysis_results" \
  --pixel-size-um 0.325
```

Replace `0.325` with the actual pixel size of the images in µm/pixel.

## Main command-line arguments

| Argument | Required | Default | Description |
|---|---:|---:|---|
| `--root-dir` | Yes | None | Root folder containing one subfolder per image. |
| `--output-dir` | Yes | None | Destination folder for all results. |
| `--pixel-size-um` | Yes | None | Pixel size in µm/pixel. Used when coordinates or geometries are in pixels. |
| `--step-um` | No | `50` | Width of each concentric ring in µm. |
| `--geojson-units` | No | `px` | Units of annotation geometries: `px` or `um`. |
| `--cells-units` | No | `um` | Units of centroid coordinates: `um` or `px`. |
| `--tumor-key` | No | `Tumor` | Keyword used to find the tumor geometry file. |
| `--ipsi-key` | No | `Cortex_ipsilateral` | Keyword used to find the ipsilateral cortex file. |
| `--contra-key` | No | `Cortex_contralateral` | Keyword used to find the contralateral cortex file. |
| `--swap-xy` | No | Off | Swaps X and Y when coordinate axes are exchanged. |
| `--flip-y` | No | Off | Flips the Y axis when exported coordinate systems use opposite Y directions. |

## Example with 50 µm rings

If GeoJSON geometries are in pixels and QuPath centroid coordinates are in µm:

```bash
python analyze_distance_normalized_ring_fast.py \
  --root-dir "/data/gephyrin_vgat/input" \
  --output-dir "/data/gephyrin_vgat/results" \
  --step-um 50 \
  --pixel-size-um 0.325 \
  --geojson-units px \
  --cells-units um
```

If both geometries and centroid coordinates are already expressed in µm:

```bash
python analyze_distance_normalized_ring_fast.py \
  --root-dir "/data/gephyrin_vgat/input" \
  --output-dir "/data/gephyrin_vgat/results" \
  --step-um 50 \
  --pixel-size-um 0.325 \
  --geojson-units um \
  --cells-units um
```

## What the Python script does

For each image and cortex side, the script:

1. Loads the tumor and cortex geometries.
2. Loads Gephyrin and vGAT centroid coordinates.
3. Keeps ipsilateral and contralateral points separate.
4. Filters points geometrically so that only detections inside the corresponding cortex are analyzed.
5. Creates concentric rings around the tumor boundary.
6. Clips each ring to the relevant cortex annotation.
7. Calculates the effective analyzed area of each clipped ring.
8. Counts Gephyrin and vGAT detections in each ring.
9. Calculates density as:

   ```text
   density_per_mm2 = count / area_um2 × 1,000,000
   ```

10. Exports per-image tables, plots, overlays, ring geometries, and PDF reports.
11. Combines all analyzed images into global summary tables.

The maximum radius is calculated automatically so that the generated buffers cover the complete cortex geometry.

---

# Output structure

Example:

```text
analysis_results/
├── Sample_01/
│   ├── Sample_01__tumor_rings_report.pdf
│   ├── ipsilateral/
│   │   ├── Sample_01__ipsilateral__rings_long.csv
│   │   ├── Sample_01__ipsilateral__rings_wide.csv
│   │   ├── Sample_01__ipsilateral__density.png
│   │   ├── Sample_01__ipsilateral__count_area.png
│   │   ├── Sample_01__ipsilateral__overlay_Gephyrin.png
│   │   ├── Sample_01__ipsilateral__overlay_vGat.png
│   │   └── Sample_01__ipsilateral__rings_um.geojson
│   └── contralateral/
│       └── ...
├── Sample_02/
│   └── ...
├── summary_rings_long_all_images.csv
├── summary_density_by_ring_mean_sem.csv
├── summary_ipsi_vs_contra_mwu_by_ring.csv
├── pooled_ipsilateral__density_mean_sem.png
└── pooled_contralateral__density_mean_sem.png
```

## Main result columns

The long-format ring table includes:

| Column | Meaning |
|---|---|
| `image` | Image folder name. |
| `cortex_side` | `ipsilateral` or `contralateral`. |
| `marker` | `Gephyrin` or `vGat`. |
| `ring_index` | Ring number starting at 0. |
| `r_in_um` | Inner ring radius in µm. |
| `r_out_um` | Outer ring radius in µm. |
| `r_mid_um` | Ring midpoint in µm. |
| `area_um2` | Ring area after clipping to the cortex. |
| `count` | Number of detections in the clipped ring. |
| `density_per_um2` | Count divided by ring area in µm². |
| `density_per_mm2` | Count normalized to one mm². |
| `cum_count` | Cumulative count from the tumor outward. |
| `cum_area_um2` | Cumulative analyzed area. |
| `cum_density_per_mm2` | Cumulative area-normalized density. |

---

# Quality control

Before interpreting the pooled results, verify the following for several representative images:

1. **Cellpose detections** correctly identify Gephyrin or vGAT puncta.
2. **Cortex classes** are correctly assigned as ipsilateral or contralateral.
3. **Tumor and cortex GeoJSON files** align spatially.
4. **Overlay PNGs** show detections inside the correct cortex.
5. **Ring contours** expand outward from the tumor and are clipped to the cortex.
6. The selected `--pixel-size-um` matches the image calibration in QuPath.
7. The `--geojson-units` and `--cells-units` arguments match the exported data.

If the geometries appear rotated or transposed, test `--swap-xy`. If they appear vertically mirrored, test `--flip-y`. These options should only be enabled after checking the overlay images.

# Common problems

## No cortex annotation selected in QuPath

Check that the annotation has the correct PathClass:

```text
Cortex_ipsilateral
```

or:

```text
Cortex_contralateral
```

## Python reports that no tumor annotation was found

Ensure the image folder contains a geometry file whose name includes the keyword `Tumor`, or provide a different keyword:

```bash
--tumor-key "MyTumorClass"
```

## Detections from the wrong side appear in an overlay

Confirm that `ipsilateral` or `contralateral` is present in each single-cell CSV filename. The script separates points by both marker and side and then applies an additional geometric cortex filter.

## Rings or detections do not align

The most common cause is a unit mismatch. Check:

- actual pixel size in QuPath;
- whether GeoJSON coordinates are in pixels or µm;
- whether CSV centroid coordinates are in pixels or µm;
- whether `--swap-xy` or `--flip-y` is required.

## An expected summary file is absent

`summary_ipsi_vs_contra_mwu_by_ring.csv` is only produced when SciPy is available and each comparison has sufficient finite values. The script requires at least three ipsilateral and three contralateral observations for a ring-level Mann-Whitney U test.

# Reproducibility recommendations

Record the following with each analysis:

- QuPath version.
- BIOP Cellpose extension version.
- Cellpose environment and model.
- Image pixel size in µm/pixel.
- Channel assignment for each marker.
- Cellpose diameter, flow threshold, and cell probability threshold.
- Ring width in µm.
- Coordinate units used for GeoJSON and centroid CSV files.
- Any use of `--swap-xy` or `--flip-y`.

