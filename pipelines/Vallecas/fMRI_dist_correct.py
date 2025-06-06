#!/usr/bin/env python3

import os, shutil, subprocess, tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn.image import smooth_img
from scipy.ndimage import binary_closing, binary_fill_holes, gaussian_filter

# ───────────────────────────── CONFIG ──────────────────────────────
ROOT_DIR = Path("bids_test")            
OUT_DIR  = "distcorr"                   
N_JOBS   = "8"                          # threads for syn ANTs
BET_FRAC = "0.3"                        # BET aggressiveness
MASK_SMOOTH_SIGMA = 3                   # σ for mask smoothing
SMOOTH_INV_MM     = 2                   # FWHM on inverted T1

# Set binary paths, this should generally not be necessary but I had some issues
os.environ["PATH"] = "/ants/install/bin:" + os.environ.get("PATH", "")
ANTSREG   = "/bin/antsRegistrationSyN.sh"
ANTSAPPLY = "/bin/antsApplyTransforms"
ANTSN4    = "/bin/N4BiasFieldCorrection"
ANTSRS    = "/bin/ResampleImageBySpacing"
ANTSIM    = "/bin/ImageMath"
BET       = "/fsl/share/fsl/bin/bet"
MCFLIRT   = "/fsl/share/fsl/bin/mcflirt"
# ───────────────────────────────────────────────────────────────────

# bin helpers ------------------------------------------------------------

def _which(cmd: str, override: str | None) -> str:
  path = override or shutil.which(cmd)
if path is None:
  raise RuntimeError(f"{cmd} not found")
return path

REG_EXE   = _which("antsRegistrationSyN.sh", ANTSREG)
APPLY_EXE = _which("antsApplyTransforms",   ANTSAPPLY)
N4_EXE    = _which("N4BiasFieldCorrection", ANTSN4) 
RS_EXE    = _which("ResampleImageBySpacing", ANTSRS)
IM_EXE    = _which("ImageMath", ANTSIM)
BET_EXE   = _which("bet", BET)
MF_EXE    = _which("bet", MCFLIRT)

# processing --------------------------------------------------------------

def make_brain_mask(
  t1: Path, 
  out_mask: Path,
  morphology: bool = True,
  struct_size: int = 15,
  pad_slices: int = 0):
  """BET + morphology"""
if out_mask.exists():
  print("mask already exists")
return out_mask
tmp = out_mask.parent / "_bet_tmp"
print("Make brain mask with FSL bet")
subprocess.run([BET_EXE, t1, tmp, "-m", "-f", BET_FRAC], check=True)
bet_mask = tmp.with_name(tmp.name + "_mask.nii.gz")

if not morphology:
  final = nib.load(bet_mask).get_fdata().astype(np.float32)

if morphology:
  print("Brain mask morphology")
gm = nib.load(bet_mask).get_fdata().astype(np.float32)
nx, ny, nz = gm.shape

if pad_slices > 0:
  gm[:, :, nz - pad_slices : nz] = 1

filled = binary_closing(gm, structure=np.ones((struct_size, struct_size, struct_size))).astype(np.float32)
filled = binary_fill_holes(filled).astype(np.float32)
smoothed = gaussian_filter(filled, sigma=MASK_SMOOTH_SIGMA)
final = (smoothed > 0.4).astype(np.float32)

nib.save(nib.Nifti1Image(final, nib.load(t1).affine), out_mask)
print("cleanup BET temp")
bet_mask.unlink(missing_ok=True)
tmp.with_suffix(".nii.gz").unlink(missing_ok=True)
return out_mask


def mask_anat(
  anat: Path, 
  mask: Path, 
  out: Path,
  inverse: bool = True,
  inv_cap: int = 2):
  
  if out.exists():
  print(f"{out} already exists")
return out
print(f"Mask and invert ({inverse}) anat image")
img = nib.load(anat)
data = img.get_fdata().astype(np.float32)
m = nib.load(mask).get_fdata() > 0
data[~m] = 0
m_img = data
if inverse:
  inv = np.zeros_like(data, np.float32)
inv[m] = 1.0 / (data[m] + np.finfo(np.float32).eps)
inv[inv > inv_cap] = 0
m_img = inv

img_masked = nib.Nifti1Image(m_img, img.affine, img.header)
img_masked = smooth_img(img_masked, SMOOTH_INV_MM)
nib.save(img_masked, out)
return out


def head_motion_corr(bold4d: Path, out_dir: Path):
  base = bold4d.name.replace("_bold.nii.gz", "")  
out_hmc = out_dir / f"{base}_desc-hmc_bold.nii.gz"
if out_hmc.exists():
  print(f"{out_hmc} already exists")
return out_hmc
else:
  subprocess.run([
    MF_EXE, 
    "-in", str(bold4d), 
    "-o", str(out_hmc), 
    "-plots", "-report"
  ], check=True)  
return out_hmc


def extract_epi_ref(
  bold4d: Path,
  out_ref: Path,
  method: str = "mean", # "mean" or "frame"
  frame_idx: int = 0,
  skip_first: int = 10,
  skip_last: int = 10):
  
  # If the reference already exists, skip re-computation
  if out_ref.exists():
  print(f"epi ref already exists: {out_ref}")
return out_ref

print(f"Extracting EPI reference ({method})")

epi_img = nib.load(str(bold4d))
data4d = epi_img.get_fdata()
affine = epi_img.affine
header = epi_img.header

n_vols = data4d.shape[-1]

if method == "mean":
  start = skip_first
end = n_vols - skip_last
if end <= start:
  print("Computing mean over all frames due to too few frames")
start = 0
end = n_vols
mean_data = np.mean(data4d[..., start:end], axis=-1)
ref_data = mean_data.astype(np.float32)

elif method == "frame":
  if frame_idx < 0 or frame_idx >= n_vols:
  raise IndexError(
    f"frame_idx={frame_idx} is out of bounds for a 4-D image with {n_vols} volumes"
  )
ref_data = data4d[..., frame_idx].astype(np.float32)

else:
  raise ValueError(f"Unknown method '{method}'. Use 'frame' or 'mean'.")
ref_img = nib.Nifti1Image(ref_data, affine, header)
nib.save(ref_img, str(out_ref))
print(f"Saved EPI reference to {out_ref}")
return out_ref

def downsample_t1_to_epi(t1_inv: Path, epi_ref: Path, out_t1_rs: Path):
  
  epi_hdr = nib.load(str(epi_ref)).header
zx, zy, zz = epi_hdr.get_zooms()[:3]  
sx, sy, sz = str(zx), str(zy), str(zz)

if not out_t1_rs.exists():
  print("Resampling T1 inverse to EPI resolution")
subprocess.run([
  RS_EXE, "3",
  str(t1_inv),
  str(out_t1_rs),
  sx, sy, sz
], check=True)

return out_t1_rs


def ants_rigid(mov: Path, fix: Path, pref: Path):
  if (pref.parent / (pref.name + "0GenericAffine.mat")).exists():
  print("affine mat already exists")
return (pref.parent / (pref.name + "0GenericAffine.mat"))
print("aligning ", mov, " to ", fix)
subprocess.run([REG_EXE, "-d", "3", "-f", fix, "-m", mov,
                "-o", pref, "-t", "r", "-p", "f", "-n", N_JOBS], check=True)


def zero_nonoverlap_slices(warped_fix: Path, epi_ref: Path):
  epi = nib.load(epi_ref).get_fdata()
fix_img = nib.load(warped_fix); fix = fix_img.get_fdata()
for z in range(epi.shape[2]):
  if np.max(epi[:, :, z]) == 0:
  fix[:, :, z] = 0
nib.save(nib.Nifti1Image(fix, fix_img.affine, fix_img.header), warped_fix)


def ants_syn(mov: Path, fix: Path, pref: Path):
  if (pref.parent / (pref.name + "1Warp.nii.gz")).exists():
  print("1Warp already exists")
return (pref.parent / (pref.name + "1Warp.nii.gz"))
print("aligning ", mov, " to ", fix)
subprocess.run([REG_EXE, "-d", "3", "-f", fix, "-m", mov,
                "-o", pref, "-t", "so", "-p", "f", "-n", N_JOBS], check=True)


def n4_bias(in_vol: Path, out_corr: Path, out_field: Path):
  if out_field.exists():
  print("bias field already exists")
return out_field
print("calculating bias field")
subprocess.run([N4_EXE, "-d", "3", "-i", in_vol,
                "-o", f"[{out_corr},{out_field}]"], check=True)


def warp_full_epi(bold4d: Path, ref_img: Path, warp_field: Path, affine_mat: Path, bias_field: Path, out_dir: Path):
  
  base = bold4d.name.replace("_bold.nii.gz", "")  
dc_4d = out_dir / f"{base}DistCorr_bold.nii.gz" 

if dc_4d.exists():
  print(dc_4d, " already exists")
return dc_4d
else:
  # register warped epi ref to bias_corrected_b0 → get bk0GenericAffine.mat
  bc_b0 = bias_field.parent / "bias_corrected_b0.nii.gz"
bk_pref = out_dir / "bk"
subprocess.run([
  REG_EXE, "-d", "3",
  "-f", str(bc_b0),
  "-m", str(ref_img),
  "-o", str(bk_pref),
  "-t", "r", "-p", "f", "-n", N_JOBS
], check=True)

bk_aff = bk_pref.with_name(bk_pref.name + "0GenericAffine.mat")

print("Warp the entire 4-D BOLD in one go")

subprocess.run([
  APPLY_EXE, "-d", "3", "-e", "3",
  "-i", str(bold4d),
  "-r", str(ref_img),
  "-o", str(dc_4d),
  "-t", str(warp_field),
  "-t", str(affine_mat),
  "-t", str(bk_aff),
  "--float",
  "-v"
],  check=True)

return dc_4d

def bias_corr_epi(dc_4d: Path, bias_field: Path, out_dir: Path):
  
  base = dc_4d.name.replace("DistCorr_bold.nii.gz", "")  
bc_path = out_dir / f"{base}DistBiasCorr_bold.nii.gz"

if bc_path.exists():
  print(bc_path, " already exists")
else:
  bias_field_img = nib.load(str(bias_field))
dc_epi = nib.load(str(dc_4d))
bias_field_data = bias_field_img.get_fdata()
data4d = dc_epi.get_fdata()
bc_epi = data4d / bias_field_data[..., None]
nib.save(nib.Nifti1Image(bc_epi, dc_epi.affine, dc_epi.header), bc_path)

return bc_path



def final_masking(bc4d: Path, out_dir: Path):
  
  base = bc4d.name.replace("_bold.nii.gz", "")  
out = out_dir / f"{base}Masked_bold.nii.gz"
subprocess.run([
  BET_EXE, 
  bc4d, 
  out, 
  # "-m", # generate binary mask
  "-f", "0.3", 
  "-F"], check=True)


# runners -----------------------------------------------------------

def process_session(
  ses: Path,
  anat_ref: str = "t2star",
  inv_anat: bool = False):
  anat, func = ses / "anat", ses / "func"

if anat_ref == "t1":
  try:
  anat_img = next(anat.glob("*_T1w.nii.gz"))
except StopIteration:
  print(f"No T1w in {ses}"); return
if anat_ref == "t2star":
  try:
  anat_img = next(anat.glob("*_T2starw.nii.gz"))
except StopIteration:
  print(f"No T1w in {ses}"); return
bolds = list(func.glob("*_bold.nii.gz"))
if not bolds:
  print(f"No BOLD runs in {ses}"); return

print(f"===== {ses.relative_to(ROOT_DIR)} =====")
print("session process: mask extraction")
if anat_ref == "t1":
  mask   = make_brain_mask(anat_img, anat / "t1_mask.nii.gz")
print("session process: mask and inverse t1")
mskd_anat = mask_anat(anat_img, mask, anat / "t1_inverted.nii.gz", inverse = inv_anat)
if anat_ref == "t2star":
  mask   = make_brain_mask(anat_img, anat / "t2star_mask.nii.gz", struct_size = 7)
print("session process: mask t2star")
mskd_anat = mask_anat(anat_img, mask, anat / "t2star_masked.nii.gz", inverse = inv_anat)


for bold in bolds:
  print("Processing", bold.name)
out_dir = bold.parent / OUT_DIR
out_dir.mkdir(exist_ok=True)

# ADD MOTION CORRECTION STEP
print("Head motion correction mcflirt")
bold_hmc = head_motion_corr(bold, out_dir)

print("Get average hmc'd epi")
epi_ref = extract_epi_ref(bold_hmc, out_dir / "epi_ref.nii.gz")

#inv_t1_rs = downsample_t1_to_epi(inv_t1, epi_ref, out_dir / "t1_inverted_rs.nii.gz")

print(f"session process: rigid registration {mskd_anat.name} -> {epi_ref.name}")
rigid_pref = out_dir / anat_ref
ants_rigid(mskd_anat, epi_ref, rigid_pref)
anat_warped = rigid_pref.parent / f"{anat_ref}Warped.nii.gz"

print("session process: zero slices outside EPI FOV")
zero_nonoverlap_slices(anat_warped, epi_ref)

print(f"session process: non‑linear SyN {epi_ref.name} -> {anat_warped.name}")
syn_pref = out_dir / "epi"
ants_syn(epi_ref, anat_warped, syn_pref)

epi_warped = syn_pref.parent / "epiWarped.nii.gz"  
warp_field = syn_pref.parent / "epi1Warp.nii.gz"
affine_mat = syn_pref.parent / "epi0GenericAffine.mat"

print("session process: N4 bias on warped epi")
bias_corr_b0 = out_dir / "bias_corrected_b0.nii.gz"
bias_field   = out_dir / "bias_field.nii.gz"
n4_bias(epi_warped, bias_corr_b0, bias_field)

print("session process: warp entire 4‑D stack")
dc4d_path = warp_full_epi(bold_hmc, epi_warped, warp_field, affine_mat,
                          bias_field, out_dir)

print("divide distortion corrected image with bias field")
bc4d_path = bias_corr_epi(dc4d_path, bias_field, out_dir)

print("session process: final brain extraction")
final_masking(bc4d_path, out_dir)


def process_dataset(root: Path):
  for subj in sorted(root.glob("sub-*")):
  ses_dirs = sorted(subj.glob("ses-*")) 
for ses in ses_dirs:
  process_session(ses, anat_ref = "t1", inv_anat = False)

if __name__ == "__main__":
  process_dataset(ROOT_DIR)

