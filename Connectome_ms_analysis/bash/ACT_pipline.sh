#!/usr/bin/env bash
set -euo pipefail

# ========================
# Config (edit if needed)
# ========================
# Root folder with subjects: INsIDER_* subdirs inside.
# Override by exporting ROOT=/path/to/INsIDER_Subj/INsIDER_Subj before running.
ROOT="${ROOT:-/home/one/Desktop/INsIDER_Subj/INsIDER_Subj}"

# Number of streamlines to generate (per subject)
STREAMLINES="${STREAMLINES:-1000000}"

# Assignment radius (mm) for robust endpoint-to-parcel matching
ASSIGN_RAD_MM="${ASSIGN_RAD_MM:-2.5}"

# Threads (auto-detect by default)
THREADS="$(nproc)"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$THREADS}"
export MRTRIX_NTHREADS="${MRTRIX_NTHREADS:-$THREADS}"

shopt -s nullglob

for SUBJ_DIR in "$ROOT"/INsIDER_*; do
  [[ -d "$SUBJ_DIR" ]] || continue
  SUBJ_BASENAME="$(basename "$SUBJ_DIR")"
  echo "== ACT :: $SUBJ_BASENAME =="

  ACTDIR="$SUBJ_DIR/ACT"
  mkdir -p "$ACTDIR"

  # Inputs
  DWI_NII="$SUBJ_DIR/dwi.nii.gz"
  BVEC="$SUBJ_DIR/bvecs.txt"
  BVAL="$SUBJ_DIR/bvals.txt"
  T1="$SUBJ_DIR/T1.nii.gz"
  FS_APARC="$SUBJ_DIR/aparc+aseg.nii.gz"
  FS_LUT="$SUBJ_DIR/FreeSurferColorLUT.txt"

  # Outputs
  DWI_MIF="$SUBJ_DIR/dwi.mif"
  DWI_MASK="$SUBJ_DIR/dwi_mask.mif"
  RESP_WM="$ACTDIR/wm_response.txt"
  RESP_GM="$ACTDIR/gm_response.txt"
  RESP_CSF="$ACTDIR/csf_response.txt"
  FOD_WM="$ACTDIR/wmfod.mif"
  FOD_GM="$ACTDIR/gm.mif"
  FOD_CSF="$ACTDIR/csf.mif"
  FIVE_TT="$ACTDIR/5tt.mif"
  GMWMI="$ACTDIR/gmwmi.mif"
  TCK="$ACTDIR/ACT_${STREAMLINES}.tck"
  # If you later want weights: SIFT2="$ACTDIR/sift2_weights.txt"
  B0_3D="$ACTDIR/b0_3d.mif"
  T1_TO_DWI_TXT="$ACTDIR/T1_to_DWI.txt"
  APARC_DWI="$ACTDIR/aparc_in_dwi.mif"
  APARC_DWI_CLAMP="$ACTDIR/aparc_in_dwi_clamped.mif"
  CONNECTOME_ACT_SUB="$ACTDIR/connectome_ACT.csv"
  CONNECTOME_ACT_ROOT="$SUBJ_DIR/connectome_ACT.csv"

  # Per-subject log
  LOG="$ACTDIR/run.log"
  exec > >(tee -a "$LOG") 2>&1

  # Sanity check
  for f in "$DWI_NII" "$BVEC" "$BVAL" "$T1" "$FS_APARC" "$FS_LUT"; do
    [[ -s "$f" ]] || { echo "  [ERROR] Missing: $f"; continue 2; }
  done

  [[ -s "$DWI_MIF" ]] || mrconvert "$DWI_NII" "$DWI_MIF" -fslgrad "$BVEC" "$BVAL" -datatype float32
  [[ -s "$DWI_MASK" ]] || dwi2mask "$DWI_MIF" "$DWI_MASK"

  if [[ ! -s "$FIVE_TT" ]]; then
    mrconvert "$FS_APARC" "$ACTDIR/aparc+aseg.mif" -force
    5ttgen freesurfer "$ACTDIR/aparc+aseg.mif" "$FIVE_TT" -lut "$FS_LUT" -force
    5ttcheck "$FIVE_TT" || true
  fi

  [[ -s "$GMWMI" ]] || 5tt2gmwmi "$FIVE_TT" "$GMWMI"

  if [[ ! -s "$RESP_WM" || ! -s "$RESP_GM" || ! -s "$RESP_CSF" ]]; then
    dwi2response dhollander "$DWI_MIF" "$RESP_WM" "$RESP_GM" "$RESP_CSF" -mask "$DWI_MASK"
  fi

  if [[ ! -s "$FOD_WM" || ! -s "$FOD_GM" || ! -s "$FOD_CSF" ]]; then
    dwi2fod msmt_csd "$DWI_MIF" \
      "$RESP_WM" "$FOD_WM" \
      "$RESP_GM" "$FOD_GM" \
      "$RESP_CSF" "$FOD_CSF" \
      -mask "$DWI_MASK"
  fi

  # Tractography (GMWMI seeding for ACT)
  if [[ ! -s "$TCK" ]]; then
    tckgen "$FOD_WM" "$TCK" \
      -act "$FIVE_TT" -backtrack -crop_at_gmwmi \
      -seed_gmwmi "$GMWMI" \
      -maxlength 250 -select "$STREAMLINES" -cutoff 0.06
  fi

  # If you want weighted connectomes later, uncomment SIFT2 + -tck_weights_in:
  # [[ -s "$SIFT2" ]] || tcksift2 "$TCK" "$FOD_WM" "$SIFT2"

  # Registration for parcellation → DWI
  if [[ ! -s "$B0_3D" ]]; then
    dwiextract "$DWI_MIF" -bzero "$ACTDIR/b0_4d.mif"
    mrconvert "$ACTDIR/b0_4d.mif" -coord 3 0 "$ACTDIR/b0_single.mif" -force
    mrconvert "$ACTDIR/b0_single.mif" -axes 0,1,2 "$B0_3D"
  fi
  [[ -s "$T1_TO_DWI_TXT" ]] || mrregister "$T1" "$B0_3D" -type rigid -rigid "$T1_TO_DWI_TXT"

  [[ -s "$APARC_DWI" ]] || mrtransform "$FS_APARC" -linear "$T1_TO_DWI_TXT" -template "$DWI_MIF" "$APARC_DWI"
  if [[ ! -s "$APARC_DWI_CLAMP" ]]; then
    mrconvert "$APARC_DWI" "$ACTDIR/aparc_in_dwi_int.mif" -datatype int32 -force
    mrcalc "$ACTDIR/aparc_in_dwi_int.mif" 0 -max "$APARC_DWI_CLAMP" -datatype int32
  fi

  # Connectome (COUNTS; fair vs TREKKER). Add -symmetric & -assignment_radial_search.
  if [[ ! -s "$CONNECTOME_ACT_SUB" ]]; then
    tck2connectome "$TCK" "$APARC_DWI_CLAMP" "$CONNECTOME_ACT_SUB" \
      -assignment_radial_search "$ASSIGN_RAD_MM" \
      -zero_diagonal -symmetric
  fi

  cp -f "$CONNECTOME_ACT_SUB" "$CONNECTOME_ACT_ROOT"
  echo "  -> ACT connectome: $CONNECTOME_ACT_ROOT"
done

echo "All ACT subjects done."
