#!/usr/bin/env bash
set -euo pipefail

# ========================
# Config (edit if needed)
# ========================
ROOT="${ROOT:-/home/one/Desktop/INsIDER_Subj/INsIDER_Subj}"   # fixed stray dot
STREAMLINES="${STREAMLINES:-1000000}"
ASSIGN_RAD_MM="${ASSIGN_RAD_MM:-2.5}"
THREADS="$(nproc)"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$THREADS}"
export MRTRIX_NTHREADS="${MRTRIX_NTHREADS:-$THREADS}"

# Helper: pick TREKKER binary per subject or global PATH
pick_trekker_bin() {
  local subj_dir="$1"
  local cand1="$subj_dir/trekker_v1.0.0-rc5"
  local cand2="$ROOT/trekker_v1.0.0-rc5"
  local cand3
  cand3="$(command -v Trekker || command -v trekker || true)"

  if [[ -f "$cand1" ]]; then chmod +x "$cand1" || true; echo "$cand1"; return; fi
  if [[ -f "$cand2" ]]; then chmod +x "$cand2" || true; echo "$cand2"; return; fi
  if [[ -n "${cand3:-}" ]]; then echo "$cand3"; return; fi

  echo "ERROR: Trekker executable not found (looked for 'trekker_v1.0.0-rc5' and 'trekker' in PATH)" >&2
  exit 3
}

shopt -s nullglob

for SUBJ_DIR in "$ROOT"/INsIDER_*; do
  [[ -d "$SUBJ_DIR" ]] || continue
  SUBJ_BASENAME="$(basename "$SUBJ_DIR")"
  echo "== TREKKER :: $SUBJ_BASENAME =="

  TREKDIR="$SUBJ_DIR/TREKKER"
  mkdir -p "$TREKDIR"

  # Inputs in each subject folder
  DWI_NII="$SUBJ_DIR/dwi.nii.gz"
  BVEC="$SUBJ_DIR/bvecs.txt"
  BVAL="$SUBJ_DIR/bvals.txt"
  T1="$SUBJ_DIR/T1.nii.gz"
  FS_APARC="$SUBJ_DIR/aparc+aseg.nii.gz"
  FS_LUT="$SUBJ_DIR/FreeSurferColorLUT.txt"

  # Shared / regenerated
  DWI_MIF="$SUBJ_DIR/dwi.mif"
  DWI_MASK="$SUBJ_DIR/dwi_mask.mif"

  # FODs
  RESP_WM="$TREKDIR/wm_response.txt"
  RESP_GM="$TREKDIR/gm_response.txt"
  RESP_CSF="$TREKDIR/csf_response.txt"
  FOD_WM="$TREKDIR/wmfod.mif"
  FOD_GM="$TREKDIR/gm.mif"
  FOD_CSF="$TREKDIR/csf.mif"

  # 5TT & GMWMI (for consistent seeding and registration reuse)
  FS_MIF="$SUBJ_DIR/aparc+aseg.mif"
  FIVE_TT_FS="$TREKDIR/5tt_freesurfer.mif"
  GMWMI="$TREKDIR/gmwmi.mif"

  # Reuse rigid from ACT if available; else compute
  T1_TO_DWI_TXT="$SUBJ_DIR/ACT/T1_to_DWI.txt"
  B0_3D="$TREKDIR/b0_3d.mif"

  # Parcellation in DWI (for connectome)
  APARC_DWI_TREK="$TREKDIR/aparc_in_dwi_trekker.mif"
  APARC_DWI_TREK_CLAMP="$TREKDIR/aparc_in_dwi_trekker_clamped.mif"

  # TREKKER IO
  FOD_WM_NII="$TREKDIR/wmfod.nii.gz"
  GMWMI_NII="$TREKDIR/gmwmi.nii.gz"
  VTK="$TREKDIR/trekker.vtk"
  TCK="$TREKDIR/trekker_${STREAMLINES}.tck"

  CONNECTOME_TREK_SUB="$TREKDIR/connectome_trekker.csv"
  CONNECTOME_TREK_ROOT="$SUBJ_DIR/connectome_trekker.csv"

  # Per-subject log
  LOG="$TREKDIR/run.log"
  exec > >(tee -a "$LOG") 2>&1

  # Sanity checks
  for f in "$DWI_NII" "$BVEC" "$BVAL" "$T1" "$FS_APARC" "$FS_LUT"; do
    [[ -s "$f" ]] || { echo "  [ERROR] Missing: $f"; continue 2; }
  done

  # Pick TREKKER binary
  TREKKER_BIN="$(pick_trekker_bin "$SUBJ_DIR")"
  echo "  Using TREKKER: $TREKKER_BIN"

  # Convert DWI to .mif (if needed) and make mask
  [[ -s "$DWI_MIF" ]] || mrconvert "$DWI_NII" "$DWI_MIF" -fslgrad "$BVEC" "$BVAL" -datatype float32
  [[ -s "$DWI_MASK" ]] || dwi2mask "$DWI_MIF" "$DWI_MASK"

  # Responses
  if [[ ! -s "$RESP_WM" || ! -s "$RESP_GM" || ! -s "$RESP_CSF" ]]; then
    dwi2response dhollander "$DWI_MIF" "$RESP_WM" "$RESP_GM" "$RESP_CSF" -mask "$DWI_MASK"
  fi

  # FODs
  if [[ ! -s "$FOD_WM" || ! -s "$FOD_GM" || ! -s "$FOD_CSF" ]]; then
    dwi2fod msmt_csd "$DWI_MIF" \
      "$RESP_WM" "$FOD_WM" \
      "$RESP_GM" "$FOD_GM" \
      "$RESP_CSF" "$FOD_CSF" \
      -mask "$DWI_MASK"
  fi

  # 5TT from FreeSurfer & GMWMI (so seeding matches ACT policy)
  [[ -s "$FS_MIF" ]] || mrconvert "$FS_APARC" "$FS_MIF" -force
  if [[ ! -s "$FIVE_TT_FS" ]]; then
    5ttgen freesurfer "$FS_MIF" "$FIVE_TT_FS" -lut "$FS_LUT" -force
    5ttcheck "$FIVE_TT_FS" || true
  fi
  [[ -s "$GMWMI" ]] || 5tt2gmwmi "$FIVE_TT_FS" "$GMWMI"

  # Registration T1->DWI (reuse from ACT if present)
  if [[ ! -s "$T1_TO_DWI_TXT" ]]; then
    if [[ ! -s "$B0_3D" ]]; then
      dwiextract "$DWI_MIF" -bzero "$TREKDIR/b0_4d.mif"
      mrconvert "$TREKDIR/b0_4d.mif" -coord 3 0 "$TREKDIR/b0_single.mif" -force
      mrconvert "$TREKDIR/b0_single.mif" -axes 0,1,2 "$B0_3D"
    fi
    mrregister "$T1" "$B0_3D" -type rigid -rigid "$T1_TO_DWI_TXT"
  fi

  # Convert inputs for TREKKER
  [[ -s "$FOD_WM_NII" ]] || mrconvert "$FOD_WM" "$FOD_WM_NII" -datatype float32
  [[ -s "$GMWMI_NII" ]] || mrconvert "$GMWMI" "$GMWMI_NII" -datatype uint8

  # TREKKER tracking (edit flags to your preferred config)
  if [[ ! -s "$VTK" ]]; then
    "$TREKKER_BIN" track "$FOD_WM_NII" \
      --seed "$GMWMI_NII" \
      --seed_count "$STREAMLINES" \
      --minDataSupport 0.05 \
      --numberOfThreads "$THREADS" \
      --verbose info \
      --output "$VTK"
  fi

  # Convert to .tck for tck2connectome
  [[ -s "$TCK" ]] || "$TREKKER_BIN" convert "$VTK" "$TCK"

  # Parcellation → DWI
  if [[ ! -s "$APARC_DWI_TREK" ]]; then
    mrtransform "$FS_APARC" -linear "$T1_TO_DWI_TXT" -template "$DWI_MIF" "$APARC_DWI_TREK"
  fi
  if [[ ! -s "$APARC_DWI_TREK_CLAMP" ]]; then
    mrconvert "$APARC_DWI_TREK" "$TREKDIR/aparc_in_dwi_trekker_int.mif" -datatype int32 -force
    mrcalc "$TREKDIR/aparc_in_dwi_trekker_int.mif" 0 -max "$APARC_DWI_TREK_CLAMP" -datatype int32
  fi

  # Connectome (COUNTS; same policy as ACT). Use same radius & symmetric.
  if [[ ! -s "$CONNECTOME_TREK_SUB" ]]; then
    tck2connectome "$TCK" "$APARC_DWI_TREK_CLAMP" "$CONNECTOME_TREK_SUB" \
      -assignment_radial_search "$ASSIGN_RAD_MM" \
      -zero_diagonal -symmetric
  fi

  cp -f "$CONNECTOME_TREK_SUB" "$CONNECTOME_TREK_ROOT"
  echo "  -> TREKKER connectome: $CONNECTOME_TREK_ROOT"
done

echo "All TREKKER subjects done."
