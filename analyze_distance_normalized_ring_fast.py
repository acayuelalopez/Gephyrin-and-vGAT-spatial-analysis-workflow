#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Análisis por rings normalizado por área (en µm) desde anotación Tumor.

FAST + SIDE-SAFE overlays
------------------------
Partiendo de tu script analyze_distance_normalized_ring_fast.py, este archivo corrige:

1) Mezcla ipsi/contra en overlays y (potencialmente) en conteos:
   - El script original agregaba los puntos por MARCADOR sin distinguir el lado.
     Si existen CSV separados (p.ej. *_vGat_ipsilateral_single_cell.csv y *_vGat_contralateral_single_cell.csv),
     el último CSV leído (según orden alfabético) SOBREESCRIBE los puntos del marcador.
     Resultado: en contralateral puedes acabar ploteando los puntos de ipsilateral (o viceversa).

   Solución:
   - Leemos puntos por (lado, marcador). El lado se infiere por nombre de fichero
     (ipsilateral/contralateral) y, si no aparece, por la columna Parent del CSV (Annotation (Cortex_* )).
   - Para cada cortex (ipsi/contra) usamos SOLO sus puntos.
   - Además, para overlays, filtramos por pertenencia geométrica al polígono cortex (contains/covers)
     como guardarraíl extra.

2) Rings no visibles en overlay:
   - El intersection() puede devolver GeometryCollection; el dibujado anterior no lo manejaba.
   - Ahora dibujamos Polygon/MultiPolygon y GeometryCollection recursivamente.
   - Cambiamos estilo de rings a líneas azules discontinuas para que se vean.

Mantiene:
- Conversión px->µm para GeoJSON con --pixel-size-um (si --geojson-units px)
- Corrección de ejes: --swap-xy y --flip-y
- Aceleraciones: max_r exp+bin, buffers cacheados, STRtree predicate contains/covers.
- Misma CLI que tu script base.

"""

import argparse
import json
import re
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

_HAS_SHAPELY = True
try:
    from shapely.geometry import shape, Polygon, MultiPolygon, Point, mapping
    from shapely.ops import unary_union
    from shapely.affinity import affine_transform, scale as shp_scale
    from shapely.strtree import STRtree
except Exception:
    _HAS_SHAPELY = False

_HAS_SCIPY = True
try:
    from scipy.stats import mannwhitneyu
except Exception:
    _HAS_SCIPY = False


# --------------------------
# Utilidades generales
# --------------------------

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def find_column(df: pd.DataFrame, contains: str, prefer_units: str = "µm") -> str:
    cand = [c for c in df.columns if contains.lower() in c.lower()]
    if not cand:
        raise ValueError(f"No se encontró columna que contenga '{contains}'")
    if prefer_units:
        for c in cand:
            if prefer_units in c:
                return c
    return cand[0]


def detect_marker(fname: str) -> Optional[str]:
    low = fname.lower()
    if "gephyrin" in low:
        return "Gephyrin"
    if "vgat" in low:
        return "vGat"
    return None


def detect_side_from_name(fname: str) -> Optional[str]:
    low = fname.lower()
    if "ipsi" in low or "ipsilateral" in low:
        return "ipsilateral"
    if "contra" in low or "contralateral" in low:
        return "contralateral"
    return None


def detect_side_from_parent_col(df: pd.DataFrame) -> Optional[str]:
    """Intenta inferir ipsi/contra desde la columna Parent (típico de QuPath export)."""
    if "Parent" not in df.columns:
        return None
    # Ejemplos: 'Annotation (Cortex_ipsilateral)' / 'Annotation (Cortex_contralateral)'
    vals = df["Parent"].astype(str).dropna().unique().tolist()
    joined = " ".join(vals).lower()
    if "cortex_ipsilateral" in joined or "cortex\\_ipsilateral" in joined:
        return "ipsilateral"
    if "cortex_contralateral" in joined or "cortex\\_contralateral" in joined:
        return "contralateral"
    return None


def read_table_like_csv(path: Path) -> pd.DataFrame:
    with open(path, "rb") as fh:
        head = fh.read(4096).decode("utf-8", errors="ignore")
    sep = "\t" if head.count("\t") > head.count(",") else ","
    return pd.read_csv(path, sep=sep, engine="c", encoding="utf-8-sig", low_memory=False)


def load_cells_csv(path: Path) -> pd.DataFrame:
    df = read_table_like_csv(path)
    col_x = find_column(df, "Centroid X")
    col_y = find_column(df, "Centroid Y")
    out = df[[col_x, col_y]].copy()
    out.columns = ["Centroid_X", "Centroid_Y"]
    out = out[np.isfinite(out["Centroid_X"]) & np.isfinite(out["Centroid_Y"])].copy()
    return out


# --------------------------
# Lectura de geometrías
# --------------------------

_num_re = re.compile(r"[-+]?\d*\.\d+|[-+]?\d+")


def _parse_flat_geojson_txt(text: str) -> Dict:
    low = text.lower()
    if "coordinates" not in low:
        raise ValueError("No se encontró token 'coordinates' en el txt")
    start = low.index("coordinates") + len("coordinates")
    end = low.find("properties", start)
    if end == -1:
        end = len(text)
    chunk = text[start:end]
    nums = [float(x) for x in _num_re.findall(chunk)]
    if len(nums) < 6:
        raise ValueError("No hay suficientes números para un polígono")
    if len(nums) % 2 != 0:
        nums = nums[:-1]
    coords = list(zip(nums[0::2], nums[1::2]))
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    return {
        "type": "FeatureCollection",
        "features": [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [coords]}, "properties": {}}],
    }


def _load_geojson_any(path: Path) -> Dict:
    txt = path.read_text(encoding="utf-8", errors="ignore").strip()
    if txt.startswith("{") or txt.startswith("["):
        return json.loads(txt)
    return _parse_flat_geojson_txt(txt)


def _load_geojson_from_zip(zip_path: Path) -> Dict:
    with zipfile.ZipFile(zip_path, "r") as zf:
        names = zf.namelist()
        cand = [n for n in names if n.lower().endswith((".geojson", ".json"))]
        if not cand:
            raise ValueError(f"No se encontró .geojson/.json dentro de {zip_path.name}")
        cand.sort(key=lambda n: (0 if n.lower().endswith(".geojson") else 1, len(n)))
        with zf.open(cand[0]) as fh:
            data = fh.read().decode("utf-8", errors="ignore")
        return json.loads(data)


def load_annotation_geometry(path: Path):
    if not _HAS_SHAPELY:
        raise RuntimeError("Falta shapely. Instala: pip install shapely")

    if path.suffix.lower() == ".zip":
        gj = _load_geojson_from_zip(path)
    else:
        gj = _load_geojson_any(path)

    geoms = []
    if isinstance(gj, dict) and gj.get("type") == "FeatureCollection":
        for f in gj.get("features", []):
            if isinstance(f, dict) and f.get("geometry"):
                geoms.append(shape(f["geometry"]))
    elif isinstance(gj, dict) and gj.get("type") == "Feature":
        geoms.append(shape(gj.get("geometry")))
    else:
        geoms.append(shape(gj))

    if not geoms:
        raise ValueError(f"No se pudieron extraer geometrías de {path.name}")

    geom = unary_union(geoms)
    try:
        if hasattr(geom, "is_valid") and (not geom.is_valid):
            geom = geom.buffer(0)
    except Exception:
        pass
    return geom


def find_annotation_file(image_dir: Path, keyword: str) -> Optional[Path]:
    exts = (".geojson", ".json", ".zip", ".txt")
    cand = [p for p in image_dir.iterdir() if p.is_file() and p.suffix.lower() in exts and keyword.lower() in p.name.lower()]
    if not cand:
        return None
    prio = {".geojson": 0, ".zip": 1, ".json": 2, ".txt": 3}
    cand.sort(key=lambda p: prio.get(p.suffix.lower(), 9))
    return cand[0]


# --------------------------
# Transformaciones (unidad + ejes)
# --------------------------

def transform_points_xy(points: np.ndarray, *, scale_factor: float = 1.0, swap_xy: bool = False,
                        flip_y: bool = False, y_ref: Optional[float] = None) -> np.ndarray:
    if points.size == 0:
        return points
    pts = points.astype(float, copy=True)
    pts *= float(scale_factor)
    if swap_xy:
        pts = pts[:, [1, 0]]
    if flip_y:
        if y_ref is None:
            raise ValueError("flip_y requiere y_ref (ymin+ymax)")
        pts[:, 1] = float(y_ref) - pts[:, 1]
    return pts


def transform_geom_xy(geom, *, scale_factor: float = 1.0, swap_xy: bool = False,
                      flip_y: bool = False, y_ref: Optional[float] = None):
    g = geom
    if scale_factor != 1.0:
        g = shp_scale(g, xfact=float(scale_factor), yfact=float(scale_factor), origin=(0.0, 0.0))
    if swap_xy:
        g = affine_transform(g, [0.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    if flip_y:
        if y_ref is None:
            raise ValueError("flip_y requiere y_ref (ymin+ymax)")
        g = affine_transform(g, [1.0, 0.0, 0.0, -1.0, 0.0, float(y_ref)])
    return g


def compute_global_yref(geoms: List, points_by_marker: Dict[str, np.ndarray]) -> float:
    ys = []
    for g in geoms:
        if g is None or getattr(g, "is_empty", False):
            continue
        b = g.bounds
        ys.extend([b[1], b[3]])
    for pts in points_by_marker.values():
        if pts.size:
            ys.extend([float(np.nanmin(pts[:, 1])), float(np.nanmax(pts[:, 1]))])
    if not ys:
        raise ValueError("No hay Y para calcular flip")
    return float(min(ys) + max(ys))


# --------------------------
# STRtree + predicados (FIX)
# --------------------------

def build_strtree(points: List[Point]):
    try:
        return STRtree(points)
    except Exception:
        return None


def _len_result(res) -> int:
    # res puede ser array de índices (numpy) o lista/array de geometrías
    try:
        return int(len(res))
    except Exception:
        return 0


def count_points_in_geom(points_list: List[Point], geom, tree=None) -> int:
    """Cuenta puntos dentro de geom.

    Shapely>=2: predicate='contains'/'covers' (query geom = polígono; tree geoms = puntos)
    """
    if not points_list:
        return 0

    if tree is not None:
        try:
            res = tree.query(geom, predicate="contains")
            n = _len_result(res)
            if n > 0:
                return n
            res2 = tree.query(geom, predicate="covers")
            return _len_result(res2)
        except TypeError:
            pass
        except Exception:
            pass

        try:
            cand = tree.query(geom)
            if len(cand) == 0:
                return 0
            cnt = 0
            for p in cand:
                if geom.contains(p):
                    cnt += 1
            if cnt == 0:
                for p in cand:
                    if geom.covers(p):
                        cnt += 1
            return cnt
        except Exception:
            pass

    cnt = 0
    for p in points_list:
        if geom.contains(p):
            cnt += 1
    if cnt == 0:
        for p in points_list:
            if geom.covers(p):
                cnt += 1
    return cnt


def indices_points_in_polygon(tree: Optional[STRtree], poly) -> Optional[np.ndarray]:
    """Devuelve índices (o None si no se puede) de puntos del STRtree dentro del polígono."""
    if tree is None:
        return None
    try:
        res = tree.query(poly, predicate="contains")
        # En shapely>=2 suele devolver ndarray de índices
        if isinstance(res, np.ndarray):
            return res
        # Si devuelve geometrías, no tenemos índices; devolvemos None.
        return None
    except Exception:
        return None


# --------------------------
# max_r acelerado (exp + bin)
# --------------------------

def _max_r_seed_from_hull(tumor_geom, cortex_geom, step_um: float) -> float:
    hull = cortex_geom.convex_hull
    if hull.is_empty:
        return float(step_um)
    coords = list(hull.exterior.coords)
    dmax = 0.0
    for x, y in coords:
        d = Point(float(x), float(y)).distance(tumor_geom)
        if d > dmax:
            dmax = d
    return float(np.ceil(dmax / step_um) * step_um)


def compute_max_r_to_cover_cortex_fast(tumor_geom, cortex_geom, step_um: float,
                                      exp_max: int = 40, bin_iters: int = 45) -> float:
    hi = max(float(step_um), _max_r_seed_from_hull(tumor_geom, cortex_geom, step_um))

    buf = tumor_geom.buffer(hi)
    it = 0
    while (not buf.covers(cortex_geom)) and it < exp_max:
        hi *= 2.0
        buf = tumor_geom.buffer(hi)
        it += 1

    if not buf.covers(cortex_geom):
        return float(np.ceil(hi / step_um) * step_um)

    lo = 0.0
    for _ in range(int(bin_iters)):
        mid = 0.5 * (lo + hi)
        if tumor_geom.buffer(mid).covers(cortex_geom):
            hi = mid
        else:
            lo = mid

    return float(np.ceil(hi / step_um) * step_um)


# --------------------------
# Rings (cache inner buffer)
# --------------------------

def make_rings_table(tumor_geom, cortex_geom,
                     points_geom_by_marker: Dict[str, List[Point]],
                     trees_by_marker: Dict[str, Optional[STRtree]],
                     step_um: float, max_radius_um: float,
                     start_at: float = 0.0) -> pd.DataFrame:

    rows = []
    n_rings = int(np.ceil((max_radius_um - start_at) / step_um))

    inner_buf = tumor_geom.buffer(float(start_at))
    for i in range(n_rings):
        r_in = start_at + i * step_um
        r_out = start_at + (i + 1) * step_um

        outer_buf = tumor_geom.buffer(float(r_out))
        ring = outer_buf.difference(inner_buf)
        inner_buf = outer_buf

        ring_clip = ring.intersection(cortex_geom)
        area = float(ring_clip.area)

        for marker, pts_list in points_geom_by_marker.items():
            if area <= 0 or ring_clip.is_empty:
                c = 0
                dens_um2 = np.nan
            else:
                c = count_points_in_geom(pts_list, ring_clip, tree=trees_by_marker.get(marker))
                dens_um2 = c / area
            rows.append({
                "ring_index": i,
                "r_in_um": r_in,
                "r_out_um": r_out,
                "r_mid_um": 0.5 * (r_in + r_out),
                "area_um2": area,
                "marker": marker,
                "count": int(c),
                "density_per_um2": float(dens_um2) if np.isfinite(dens_um2) else np.nan,
                "density_per_mm2": float(dens_um2 * 1e6) if np.isfinite(dens_um2) else np.nan,
            })

    df = pd.DataFrame(rows)
    df["cum_count"] = df.groupby("marker")["count"].cumsum()
    df["cum_area_um2"] = df.groupby("marker")["area_um2"].cumsum()
    df["cum_density_per_um2"] = df["cum_count"] / df["cum_area_um2"].replace(0, np.nan)
    df["cum_density_per_mm2"] = df["cum_density_per_um2"] * 1e6
    return df


# --------------------------
# PLOTS
# --------------------------

def plot_density_by_ring(df_long: pd.DataFrame, title: str, out_png: Path, out_pdf: Optional[PdfPages] = None):
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    for marker, sub in df_long.groupby("marker"):
        x = sub["r_mid_um"].to_numpy()
        y = sub["density_per_mm2"].to_numpy()
        ax.plot(x, y, marker="o", ms=3, lw=1.5, label=marker)
    ax.set_xlabel("Distancia desde el borde de Tumor (µm)")
    ax.set_ylabel("Densidad (cells/mm²)")
    ax.set_title(title)
    ax.legend()
    fig.savefig(out_png, dpi=200)
    if out_pdf is not None:
        out_pdf.savefig(fig)
    plt.close(fig)


def plot_count_area(df_long: pd.DataFrame, title: str, out_png: Path, out_pdf: Optional[PdfPages] = None):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    ax = axes[0]
    for marker, sub in df_long.groupby("marker"):
        x = sub["r_mid_um"].to_numpy()
        y = sub["count"].to_numpy()
        ax.plot(x, y, marker="o", ms=3, lw=1.5, label=marker)
    ax.set_xlabel("Distancia desde el borde de Tumor (µm)")
    ax.set_ylabel("Count (cells)")
    ax.set_title("Conteo por ring")
    ax.legend()

    ax = axes[1]
    sub0 = df_long[df_long["marker"] == df_long["marker"].iloc[0]].copy()
    x = sub0["r_mid_um"].to_numpy()
    y = sub0["area_um2"].to_numpy()
    ax.plot(x, y, marker="o", ms=3, lw=1.5, color="black")
    ax.set_xlabel("Distancia desde el borde de Tumor (µm)")
    ax.set_ylabel("Área ring (µm²)")
    ax.set_title("Área por ring (recortado al cortex)")

    fig.suptitle(title)
    fig.savefig(out_png, dpi=200)
    if out_pdf is not None:
        out_pdf.savefig(fig)
    plt.close(fig)


def _draw_geom(ax, g, **kw):
    """Dibuja contorno de Polygon/MultiPolygon/GeometryCollection."""
    if g is None or getattr(g, "is_empty", False):
        return
    gt = getattr(g, "geom_type", None)
    if gt == "Polygon":
        xs, ys = g.exterior.xy
        ax.plot(np.asarray(xs), np.asarray(ys), **kw)
    elif gt == "MultiPolygon":
        for p in g.geoms:
            _draw_geom(ax, p, **kw)
    elif gt == "GeometryCollection":
        for gg in g.geoms:
            _draw_geom(ax, gg, **kw)
    else:
        # intentar acceder a .geoms si existe
        if hasattr(g, "geoms"):
            for gg in g.geoms:
                _draw_geom(ax, gg, **kw)


def plot_overlay_one_marker(tumor_geom, cortex_geom, pts: np.ndarray, marker_name: str, step_um: float,
                           out_png: Path, out_pdf: Optional[PdfPages] = None, n_rings_draw: int = 8):
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)

    _draw_geom(ax, cortex_geom, color="#aaaaaa", lw=1.2, label="Cortex")
    _draw_geom(ax, tumor_geom, color="#d62728", lw=1.6, label="Tumor")

    # rings (contornos): líneas azules discontinuas
    for i in range(n_rings_draw):
        r = (i + 1) * step_um
        ring_outer = tumor_geom.buffer(r)
        ring_clip = ring_outer.intersection(cortex_geom)
        _draw_geom(ax, ring_clip, color="black", lw=1.0, ls="--", alpha=0.9)

    colors = {"Gephyrin": "#1f77b4", "vGat": "#ff7f0e"}
    if pts is not None and pts.size:
        ax.scatter(pts[:, 0], pts[:, 1], s=1, alpha=0.30, edgecolors="none",
                   color=colors.get(marker_name, "black"), label=marker_name)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title(f"Overlay: Tumor, Cortex y rings (µm) · {marker_name}")
    from collections import OrderedDict

    handles, labels = ax.get_legend_handles_labels()
    uniq = OrderedDict()
    for h, l in zip(handles, labels):
        if not l or l.startswith("_"):
            continue
        if l not in uniq:
            uniq[l] = h

    ax.legend(uniq.values(), uniq.keys(), markerscale=2, fontsize=8)
    fig.savefig(out_png, dpi=200)
    if out_pdf is not None:
        out_pdf.savefig(fig)
    plt.close(fig)


# --------------------------
# Carga de puntos por lado+marcador
# --------------------------

def load_points_by_side_marker(image_dir: Path,
                               pixel_size_um: float,
                               cells_units: str,
                               swap_xy: bool,
                               flip_y: bool,
                               y_ref: Optional[float]) -> Tuple[
                                   Dict[str, Dict[str, np.ndarray]],
                                   Dict[str, Dict[str, List[Point]]],
                                   Dict[str, Dict[str, Optional[STRtree]]]
                               ]:
    """Lee todos los CSV single_cell y devuelve:

    - points_xy[side][marker] = Nx2 (float)
    - points_geom[side][marker] = List[Point]
    - trees[side][marker] = STRtree o None

    side ∈ {'ipsilateral','contralateral','unknown'}.

    Nota: si flip_y=True, y_ref debe pasarse ya calculado.
    """
    sides = ["ipsilateral", "contralateral", "unknown"]
    markers = ["Gephyrin", "vGat"]

    points_xy = {s: {m: np.zeros((0, 2)) for m in markers} for s in sides}
    points_geom = {s: {m: [] for m in markers} for s in sides}
    trees = {s: {m: None for m in markers} for s in sides}

    cell_csvs = sorted([p for p in image_dir.glob("*.csv") if "single_cell" in p.name.lower() and "summary" not in p.name.lower()])
    for csvp in cell_csvs:
        marker = detect_marker(csvp.name)
        if marker is None:
            continue

        df_raw = read_table_like_csv(csvp)
        side = detect_side_from_name(csvp.name) or detect_side_from_parent_col(df_raw) or "unknown"

        dfc = load_cells_csv(csvp)  # solo XY
        pts = dfc[["Centroid_X", "Centroid_Y"]].to_numpy(float)

        s_cells = pixel_size_um if cells_units == "px" else 1.0
        pts = transform_points_xy(pts, scale_factor=s_cells, swap_xy=swap_xy, flip_y=False)
        if flip_y:
            pts = transform_points_xy(pts, flip_y=True, y_ref=y_ref)

        # concatenar (por si hay varios CSV del mismo lado/marcador)
        if points_xy[side][marker].size == 0:
            points_xy[side][marker] = pts
        else:
            points_xy[side][marker] = np.vstack([points_xy[side][marker], pts])

    # construir geom + trees
    for s in sides:
        for m in markers:
            pts = points_xy[s][m]
            pts_list = [Point(float(x), float(y)) for x, y in pts] if pts.size else []
            points_geom[s][m] = pts_list
            trees[s][m] = build_strtree(pts_list) if pts_list else None

    return points_xy, points_geom, trees


# --------------------------
# Procesado por imagen
# --------------------------

def process_image_dir(image_dir: Path, out_dir: Path, step_um: float, keywords: Dict[str, str], *,
                      pixel_size_um: float,
                      geojson_units: str,
                      cells_units: str,
                      swap_xy: bool,
                      flip_y: bool) -> List[pd.DataFrame]:

    if not _HAS_SHAPELY:
        raise RuntimeError("Falta shapely. Instala: pip install shapely")

    image_name = image_dir.name
    out_img = ensure_dir(out_dir / image_name)

    # Anotaciones (para y_ref y para análisis)
    tumor_file = find_annotation_file(image_dir, keywords["tumor"])
    ipsi_file = find_annotation_file(image_dir, keywords["ipsi"])
    contra_file = find_annotation_file(image_dir, keywords["contra"])

    if tumor_file is None:
        print(f"[WARN] {image_name}: no encontré anotación Tumor", file=sys.stderr)
        return []

    tumor_raw = load_annotation_geometry(tumor_file)
    s_geo = pixel_size_um if geojson_units == "px" else 1.0
    tumor_geom = transform_geom_xy(tumor_raw, scale_factor=s_geo, swap_xy=swap_xy, flip_y=False)

    cortex_geoms = {}
    for side, cf in (("ipsilateral", ipsi_file), ("contralateral", contra_file)):
        if cf is None:
            continue
        g_raw = load_annotation_geometry(cf)
        cortex_geoms[side] = transform_geom_xy(g_raw, scale_factor=s_geo, swap_xy=swap_xy, flip_y=False)

    # y_ref: lo calculamos con geometrías (sin puntos todavía) usando su bbox
    if flip_y:
        ys = []
        for g in [tumor_geom] + list(cortex_geoms.values()):
            if g is None or getattr(g, 'is_empty', False):
                continue
            b = g.bounds
            ys.extend([b[1], b[3]])
        if not ys:
            raise ValueError('No hay Y para calcular flip')
        y_ref = float(min(ys) + max(ys))
        tumor_geom = transform_geom_xy(tumor_geom, flip_y=True, y_ref=y_ref)
        for side in list(cortex_geoms.keys()):
            cortex_geoms[side] = transform_geom_xy(cortex_geoms[side], flip_y=True, y_ref=y_ref)
    else:
        y_ref = None

    # Cargar puntos por lado+marcador (ya con flip si aplica)
    points_xy, points_geom, trees = load_points_by_side_marker(
        image_dir,
        pixel_size_um=pixel_size_um,
        cells_units=cells_units,
        swap_xy=swap_xy,
        flip_y=flip_y,
        y_ref=y_ref,
    )

    dfs_out = []
    report_pdf_path = out_img / f"{image_name}__tumor_rings_report.pdf"

    with PdfPages(report_pdf_path) as pdf:
        for side, cortex_geom in cortex_geoms.items():
            max_r = compute_max_r_to_cover_cortex_fast(tumor_geom, cortex_geom, step_um=step_um)

            # Seleccionar puntos del lado (si no hay, usar unknown filtrado por geometría)
            pts_xy_side = {m: points_xy.get(side, {}).get(m, np.zeros((0, 2))) for m in ("Gephyrin", "vGat")}
            pts_geom_side = {m: points_geom.get(side, {}).get(m, []) for m in ("Gephyrin", "vGat")}
            trees_side = {m: trees.get(side, {}).get(m, None) for m in ("Gephyrin", "vGat")}

            # fallback si no hay puntos explícitos de ese lado
            if all(pts_xy_side[m].size == 0 for m in pts_xy_side):
                pts_xy_side = {m: points_xy["unknown"][m] for m in ("Gephyrin", "vGat")}
                pts_geom_side = {m: points_geom["unknown"][m] for m in ("Gephyrin", "vGat")}
                trees_side = {m: trees["unknown"][m] for m in ("Gephyrin", "vGat")}

            # Guardarraíl: filtrar por cortex_geom para overlays (evita mezclas visuales)
            pts_xy_side_filtered = {}
            pts_geom_side_filtered = {}
            trees_side_filtered = {}
            for m in ("Gephyrin", "vGat"):
                pts = pts_xy_side[m]
                pts_list = pts_geom_side[m]
                tree = trees_side[m]
                if pts.size and tree is not None:
                    idx = indices_points_in_polygon(tree, cortex_geom)
                    if idx is not None:
                        pts_xy_side_filtered[m] = pts[idx]
                        pts_geom_side_filtered[m] = [pts_list[int(i)] for i in idx]
                        trees_side_filtered[m] = build_strtree(pts_geom_side_filtered[m]) if len(idx) else None
                    else:
                        # fallback: no índices -> filtrar a mano (más lento pero solo para overlay)
                        mask = np.array([cortex_geom.covers(p) for p in pts_list], dtype=bool)
                        pts_xy_side_filtered[m] = pts[mask]
                        pts_geom_side_filtered[m] = [p for p, keep in zip(pts_list, mask) if keep]
                        trees_side_filtered[m] = build_strtree(pts_geom_side_filtered[m]) if pts_geom_side_filtered[m] else None
                else:
                    pts_xy_side_filtered[m] = pts
                    pts_geom_side_filtered[m] = pts_list
                    trees_side_filtered[m] = tree

            # Cuantificación usando puntos filtrados al lado (más consistente con lo que quieres)
            df_long = make_rings_table(
                tumor_geom,
                cortex_geom,
                pts_geom_side_filtered,
                trees_side_filtered,
                step_um=step_um,
                max_radius_um=max_r,
            )

            df_long.insert(0, "image", image_name)
            df_long.insert(1, "cortex_side", side)
            df_long.insert(2, "tumor_file", tumor_file.name)
            df_long.insert(3, "cortex_file", (ipsi_file.name if side == "ipsilateral" and ipsi_file else (contra_file.name if contra_file else "")))
            df_long.insert(4, "pixel_size_um", float(pixel_size_um))
            df_long.insert(5, "geojson_units", geojson_units)
            df_long.insert(6, "cells_units", cells_units)
            df_long.insert(7, "swap_xy", bool(swap_xy))
            df_long.insert(8, "flip_y", bool(flip_y))

            out_side = ensure_dir(out_img / side)
            df_long.to_csv(out_side / f"{image_name}__{side}__rings_long.csv", index=False)

            df_wide = df_long.pivot_table(
                index=["ring_index", "r_in_um", "r_out_um", "r_mid_um", "area_um2"],
                columns="marker",
                values=["count", "density_per_mm2", "cum_count", "cum_density_per_mm2"],
            )
            df_wide.columns = [f"{a}_{b}" for a, b in df_wide.columns]
            df_wide = df_wide.reset_index()
            df_wide.to_csv(out_side / f"{image_name}__{side}__rings_wide.csv", index=False)

            plot_density_by_ring(df_long, f"{image_name} · {side} · densidad (count/area)", out_side / f"{image_name}__{side}__density.png", pdf)
            plot_count_area(df_long, f"{image_name} · {side}", out_side / f"{image_name}__{side}__count_area.png", pdf)

            # Overlays por marcador, usando puntos filtrados
            for marker_name in ("Gephyrin", "vGat"):
                plot_overlay_one_marker(
                    tumor_geom,
                    cortex_geom,
                    pts_xy_side_filtered.get(marker_name, np.zeros((0, 2))),
                    marker_name,
                    step_um=step_um,
                    out_png=out_side / f"{image_name}__{side}__overlay_{marker_name}.png",
                    out_pdf=pdf,
                    n_rings_draw=8,
                )

            # Export rings en µm
            feats = []
            for ri in sorted(df_long["ring_index"].unique()):
                r_out = float(df_long.loc[df_long["ring_index"] == ri, "r_out_um"].iloc[0])
                r_in = float(df_long.loc[df_long["ring_index"] == ri, "r_in_um"].iloc[0])
                ring = tumor_geom.buffer(r_out).difference(tumor_geom.buffer(r_in)).intersection(cortex_geom)
                if ring.is_empty:
                    continue
                feats.append({
                    "type": "Feature",
                    "properties": {"ring_index": int(ri), "r_in_um": r_in, "r_out_um": r_out, "cortex_side": side, "image": image_name, "units": "um"},
                    "geometry": mapping(ring),
                })
            if feats:
                gj = {"type": "FeatureCollection", "features": feats}
                (out_side / f"{image_name}__{side}__rings_um.geojson").write_text(json.dumps(gj), encoding="utf-8")

            dfs_out.append(df_long)

    return dfs_out


# --------------------------
# Resumen global
# --------------------------

def summarize_across_images(df_all: pd.DataFrame, out_dir: Path):
    ensure_dir(out_dir)
    df_all.to_csv(out_dir / "summary_rings_long_all_images.csv", index=False)

    grp = df_all.groupby(["cortex_side", "marker", "ring_index", "r_mid_um"], as_index=False)
    agg = grp.agg(
        n_images=("image", "nunique"),
        mean_density_mm2=("density_per_mm2", "mean"),
        sem_density_mm2=("density_per_mm2", lambda x: np.nanstd(x, ddof=1) / np.sqrt(np.sum(np.isfinite(x))) if np.sum(np.isfinite(x)) > 1 else np.nan),
        mean_count=("count", "mean"),
        mean_area_um2=("area_um2", "mean"),
    )
    agg.to_csv(out_dir / "summary_density_by_ring_mean_sem.csv", index=False)

    for side in sorted(agg["cortex_side"].unique()):
        fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
        sub_side = agg[agg["cortex_side"] == side]
        for marker in sorted(sub_side["marker"].unique()):
            sub = sub_side[sub_side["marker"] == marker].sort_values("r_mid_um")
            x = sub["r_mid_um"].to_numpy(); m = sub["mean_density_mm2"].to_numpy(); s = sub["sem_density_mm2"].to_numpy()
            ax.plot(x, m, label=marker)
            ax.fill_between(x, m - s, m + s, alpha=0.2)
        ax.set_xlabel("Distancia desde el borde de Tumor (µm)")
        ax.set_ylabel("Densidad media (cells/mm²)")
        ax.set_title(f"Pooled · {side} · densidad media ± SEM")
        ax.legend()
        fig.savefig(out_dir / f"pooled_{side}__density_mean_sem.png", dpi=200)
        plt.close(fig)

    if _HAS_SCIPY:
        try:
            rows = []
            for marker in df_all["marker"].unique():
                for ri in sorted(df_all["ring_index"].unique()):
                    a = df_all[(df_all["marker"] == marker) & (df_all["cortex_side"] == "ipsilateral") & (df_all["ring_index"] == ri)]
                    b = df_all[(df_all["marker"] == marker) & (df_all["cortex_side"] == "contralateral") & (df_all["ring_index"] == ri)]
                    if a.empty or b.empty:
                        continue
                    x = a["density_per_mm2"].to_numpy(float); y = b["density_per_mm2"].to_numpy(float)
                    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
                    if len(x) < 3 or len(y) < 3:
                        continue
                    u = mannwhitneyu(x, y, alternative="two-sided")
                    rows.append({"marker": marker, "ring_index": int(ri), "n_ipsi": int(len(x)), "n_contra": int(len(y)), "mwu_stat": float(u.statistic), "mwu_p": float(u.pvalue)})
            if rows:
                pd.DataFrame(rows).to_csv(out_dir / "summary_ipsi_vs_contra_mwu_by_ring.csv", index=False)
        except Exception:
            pass


# --------------------------
# Main
# --------------------------

def main():
    ap = argparse.ArgumentParser(description="Rings de 50 µm desde Tumor (normalizado por área) FAST + side-safe overlays.")
    ap.add_argument("--root-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--step-um", type=float, default=50.0)

    ap.add_argument("--pixel-size-um", type=float, required=True, help="µm/px")
    ap.add_argument("--geojson-units", choices=["px", "um"], default="px")
    ap.add_argument("--cells-units", choices=["um", "px"], default="um")

    ap.add_argument("--swap-xy", action="store_true")
    ap.add_argument("--flip-y", action="store_true")

    ap.add_argument("--tumor-key", default="Tumor")
    ap.add_argument("--ipsi-key", default="Cortex_ipsilateral")
    ap.add_argument("--contra-key", default="Cortex_contralateral")

    args = ap.parse_args()

    root_dir = Path(args.root_dir)
    out_dir = ensure_dir(Path(args.output_dir))

    keywords = {"tumor": args.tumor_key, "ipsi": args.ipsi_key, "contra": args.contra_key}

    all_dfs = []
    for image_dir in sorted([p for p in root_dir.iterdir() if p.is_dir()]):
        try:
            dfs = process_image_dir(
                image_dir,
                out_dir,
                step_um=float(args.step_um),
                keywords=keywords,
                pixel_size_um=float(args.pixel_size_um),
                geojson_units=args.geojson_units,
                cells_units=args.cells_units,
                swap_xy=bool(args.swap_xy),
                flip_y=bool(args.flip_y),
            )
            all_dfs.extend(dfs)
        except Exception as e:
            print(f"[WARN] {image_dir.name}: fallo procesando: {e}", file=sys.stderr)

    if all_dfs:
        df_all = pd.concat(all_dfs, ignore_index=True)
        summarize_across_images(df_all, out_dir)
        print(f"[OK] Terminado. Resultados en: {out_dir}")
    else:
        print("[WARN] No se generaron resultados (¿faltan anotaciones o CSVs?)", file=sys.stderr)


if __name__ == "__main__":
    main()