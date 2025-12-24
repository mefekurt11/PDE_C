# PDE-Based Image Compression (Edge-Based + Homogeneous Diffusion)

**Workspace:** `C_implementation`  
**Date:** 24 Aralık 2025  
**Language:** C (C11)  
**External libs:** `stb_image.h`, `stb_image_write.h`, `jbigkit` (`jbig.h`, `libjbig`)  

This document describes what was implemented in the project, which parameters were used, and the mathematical background of the method. The implementation follows the general pipeline from Mainberger et al. (edge-based image compression with homogeneous diffusion) with practical engineering adjustments that were needed to get good reconstructions on real images.

---

## 1) What we built (high-level)

We implemented an **encoder/decoder** that compresses an image by storing only:

1. A **binary mask** of “known pixels” (initially edges + image border; later extended with optional **background anchor samples**).
2. A **1‑D signal** containing the RGB values at those known pixels (in raster order).
3. A reconstruction stage that fills all unknown pixels by solving the **Laplace equation** (homogeneous diffusion) with **Dirichlet boundary conditions** imposed at the known pixels.

Reconstruction is done independently for R, G, and B channels.

---

## 2) Current repository structure (relevant files)

- `main.c`
  - CLI entry point.
  - Supports `demo`, `encode`, `decode`.
  - Defines `STB_IMAGE_IMPLEMENTATION` and `STB_IMAGE_WRITE_IMPLEMENTATION`.

- `pde_codec.h`
  - Public API and `Grid` structure.

- `pde_codec.c`
  - Core implementation:
    - Gaussian presmoothing (2D and 1D)
    - Marr-Hildreth style edge detection (practical variant)
    - JBIG encode/decode for binary masks
    - Color signal extraction/scatter
    - Max–Lloyd quantization
    - Multigrid-style solver (restriction + relaxation + W‑cycle)
    - **Background sampling / anchor points** (important quality fix)
    - File-based encode/decode (`.pde` container)

- `Makefile`
  - Builds `pde_codec` and links JBIG and libm.

---

## 3) The mathematical model

### 3.1 Missing-data view of the background

Let the image domain be a discrete grid:

- \(\Omega = \{(x,y) : 0 \le x < W,\ 0 \le y < H\}\)

Let \(K \subset \Omega\) be the set of **known pixels** (stored pixels). In the paper, \(K\) contains semantic edges and the image boundary (frame). In our implementation, \(K\) is:

- **Edges** from edge detection
- **Image border** (implicitly included by extraction rule)
- **Optional background anchors** inside large uniform background regions to avoid color drift

For each color channel \(u\) (R/G/B), reconstruction is:

\[
\begin{cases}
\Delta u(x,y) = 0, & (x,y) \in \Omega \setminus K\\
 u(x,y) = f(x,y), & (x,y) \in K
\end{cases}
\]

This is a **Dirichlet boundary value problem** for the discrete Laplacian.

Interpretation: all pixels not stored are **missing data**, and the background is obtained as the “smoothest possible” interpolation between known values.

### 3.2 Discrete Laplacian (5-point stencil)

On a regular grid, the Laplacian is approximated by:

\[
(\Delta u)_{i,j} \approx u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}
\]

The homogeneous diffusion equation \(\Delta u = 0\) means:

\[
 u_{i,j} = \frac{1}{4} \left(u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}\right)
\]

So each unknown pixel becomes the average of its neighbors.

### 3.3 Dirichlet conditions (anchoring)

Pixels in \(K\) are fixed:

\[
 u_{i,j} = f_{i,j} \quad \text{for}\ (i,j) \in K
\]

In code, this is `mask[i] == 1` ⇒ `u[i] = f[i]`.

**Why background can fail without anchors:** If the background is a large region far from edges, the reconstruction is influenced by distant boundaries and edges. If edges have strong colors, diffusion can create banding or color bleeding. Adding a few **interior Dirichlet anchors** stabilizes the solution.

---

## 4) Algorithmic pipeline (what happens in our codec)

### 4.1 Encoder pipeline (implemented in `pde_encode()`)

1. **Load input PNG** as RGB (8-bit) via stb.
2. **Split channels** into float arrays `R`,`G`,`B`.
3. **Presmooth for edge detection only**:
   - Make copies `R_smooth`,`G_smooth`,`B_smooth`.
   - Apply separable Gaussian blur `gaussian_blur_2d()`.

4. **Edge detection** with `detect_edges_marr()`:
   - Computes an approximation of the Laplacian magnitude across channels.
   - Thresholds it to produce a binary edge mask.

5. **Background anchoring (extra Dirichlet points)**:
   - `add_background_samples(...)`: grid-like sampling away from edges.
   - `detect_uniform_regions(...)`: block variance test, anchors center of uniform blocks.

6. **Encode the mask with JBIG**:
   - Mask is packed to 1-bit bitmap and compressed with JBIG.

7. **Extract 1‑D color signals** in raster order:
   - `extract_color_signal()` returns values at all `mask==1` plus border pixels.
   - Also returns `positions[]` array of indices to map signal ↔ image.

8. **Optional 1‑D presmoothing**:
   - `gaussian_blur_1d(signal, ..., presmooth_sigma * 0.5)`.

9. **Quantize signal** (Max–Lloyd):
   - `quantize_channel(signal, len, quant_levels)`.

10. **Write `.pde` file**:
   - Header (magic, W,H, sigma, threshold, q)
   - JBIG blob length + bytes
   - signal length + RGB triplets, stored as 8-bit per sample (after quantization)


### 4.2 Decoder pipeline (implemented in `pde_decode()`)

1. Read `.pde` header.
2. Decode JBIG mask.
3. Read 1‑D signal RGB and rebuild float signal arrays.
4. Reconstruct `positions[]` by scanning raster order and selecting pixels where `(mask==1) OR border`.
5. Create 3 grids (`Grid` struct) for R/G/B.
6. Put known values into `f` and `u`, and set `mask`.
7. Solve \(\Delta u=0\) on unknown pixels using multigrid-ish W-cycles + relaxation.
8. Write output PNG.

---

## 5) Edge detection details

### 5.1 Presmoothing (Gaussian)

We use separable Gaussian blur in 2D:

\[
G_\sigma(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-x^2/(2\sigma^2)}
\]

Kernel radius is chosen as \(r = \lceil 3\sigma \rceil\). Convolution is applied as horizontal then vertical 1‑D passes.

### 5.2 Practical Marr-Hildreth variant

Classic Marr-Hildreth uses Laplacian of Gaussian and looks for zero crossings. Our implementation uses a simplified Laplacian magnitude threshold. For each interior pixel:

\[
L_R = \Delta R,\ L_G = \Delta G,\ L_B = \Delta B
\]

Combine to one scalar:

\[
L = \sqrt{L_R^2 + L_G^2 + L_B^2}
\]

Edge decision:

\[
\text{edge} \iff L > T
\]

This is a pragmatic approximation (works well for cartoon/illustration style images) and is simpler than true zero-crossing logic.

---

## 6) Background as missing data + boundary anchoring (what we changed)

### 6.1 The problem we observed

On images with large uniform backgrounds (like Spider-Man on white/gray), using only “semantic edges + border” as Dirichlet points can cause:

- Background gradients that should not exist
- Color bleeding from object edges into the background

Mathematically, even though \(\Delta u=0\) is correct, the solution is entirely defined by boundary constraints. If constraints are sparse and far away, the harmonic function that fits them may still produce unwanted smooth variations.

### 6.2 The fix we implemented

We added **extra Dirichlet constraints** in the background:

1. **Regular background sampling** (`add_background_samples`)
   - For points on a coarse grid (spacing `bg_sample_spacing`), add a point if there is no edge nearby.

2. **Uniform-block sampling** (`detect_uniform_regions`)
   - Partition image into blocks.
   - Compute color variance in each block.
   - If variance < threshold and the block contains no edges, anchor the block center.

That converts “background” into a better constrained missing-data problem.

### 6.3 Quality impact (measured)

After adding background anchors and tuning default parameters, demo output achieved:

- **PSNR ≈ 33.65 dB**
- **Bitrate ≈ 4.01 bpp** on the tested image

(Values printed by the demo after changes.)

---

## 7) Quantization (Max–Lloyd)

We quantize only the extracted 1‑D signals (not the entire image).

Given a signal \(s\) and desired \(q\) levels, Max–Lloyd iteratively updates:

- decision boundaries \(b_i\)
- reconstruction points \(r_i\)

Boundaries:

\[
 b_i = \frac{r_i + r_{i+1}}{2}
\]

Assignment step: pick bin index for each sample.

Update step:

\[
 r_i \leftarrow \frac{1}{|S_i|}\sum_{x\in S_i} x
\]

Stop when reconstruction points change less than a small threshold.

In the file format, we store quantized samples as 8-bit integers.

---

## 8) Multigrid-style solver (implementation details)

### 8.1 Relaxation (Gauss-Seidel-like)

For each unknown pixel:

\[
 u_{i,j} \leftarrow \frac{1}{N}\sum_{(p,q)\in \mathcal{N}(i,j)} u_{p,q}
\]

Where \(\mathcal{N}\) is the 4-neighborhood and \(N\) is number of valid neighbors (edges have fewer neighbors).

### 8.2 Restriction (coarsening)

We build a coarse grid by grouping 2×2 blocks. For coarse pixel:

- If any fine pixel is a Dirichlet point, coarse becomes Dirichlet too and averages the known `f` values.
- Otherwise, coarse `u` is the average of fine `u`.

### 8.3 W-cycle

At each level:

1. Pre-smooth `relax(g, 3)`
2. Restrict to coarse
3. Two recursive calls to coarse (`W` shape)
4. Prolong/assign coarse to fine (current code uses nearest mapping for unknown pixels)
5. Post-smooth `relax(g, 4)`

Decoder runs multiple W-cycles:

- 3 W-cycles + final relaxation (`relax(..., 50)`) for convergence.

---

## 9) Parameters used (defaults and tuned values)

### 9.1 Default CLI params (as last set in `main.c`)

- `presmooth_sigma = 0.5`
- `edge_threshold = 8.0`
- `quant_levels = 128`
- `mg_levels = 10`

### 9.2 Background anchoring parameters (current hardcoded values)

Inside `pde_encode()`:

- `bg_sample_spacing = 16`
- `variance_threshold = 100.0`
- `block_size = 24`

These values were tuned to improve background reconstruction on uniform backgrounds.

### 9.3 Previously tested / observed values during tuning

- Presmoothing sigma: 0.8 → 0.5 (less blur improves edge localization)
- Threshold sweep (with q=256, sigma=0.5):
  - `t=5.0` produced ~27 dB PSNR but higher bpp
  - `t≈8..12` produced slightly lower PSNR but reduced bpp

---

## 10) File format (`.pde`) used

Binary layout (little-endian as written by C on macOS):

1. Magic: 4 bytes: `"PDE1"`
2. `int w`, `int h`
3. `float presmooth_sigma`, `float edge_threshold`
4. `int quant_levels`
5. `size_t jbig_len`
6. `jbig_len` bytes mask data
7. `int sig_len`
8. `sig_len` × (R,G,B) 8-bit samples

Note: `size_t` makes format platform-dependent (8 bytes on 64-bit macOS). If portability is required, change to fixed-width `uint64_t`.

---

## 11) Build configuration

`Makefile`:

- `CC = gcc`
- `CFLAGS = -O3 -std=c11 -I/opt/homebrew/include -Wall`
- `LDFLAGS = -L/opt/homebrew/lib -ljbig -lm`

Build:

- `make clean && make`

---

## 12) Key engineering decisions / fixes we made

1. **Corrected build to C** (Makefile used `.cpp`; repo uses `.c`).
2. **Avoided including stb in the public header** to prevent redefinition errors.
3. **Fixed encoder to extract colors from the original image**, not the presmoothed version.
4. **Reduced signal storage from 16-bit to 8-bit** (after quantization), significantly decreasing bpp.
5. **Improved convergence** by running multiple W-cycles + final relaxation.
6. **Major quality fix:** background anchor sampling (Dirichlet points) to stabilize large uniform areas.

---

## 13) Known limitations and next improvements

1. **Edge detection is simplified** (threshold on Laplacian magnitude). True Marr-Hildreth would do proper zero-crossing logic.
2. **Prolongation isn’t fully used** in W-cycle (current fine update is a nearest mapping). Using `prolongate()` consistently could improve convergence/quality.
3. **No entropy coding for 1-D signal** yet (paper uses predictive + entropy coding). Adding Huffman/ANS or zstd would reduce bpp.
4. **File format portability**: replace `size_t` with fixed width.
5. **Parameter exposure**: background sample spacing/variance threshold are hard-coded; should be CLI parameters.

---

## 14) How to run

From the project folder:

- Build:
  - `make clean && make`

- Demo:
  - `./pde_codec demo`

- Encode:
  - `./pde_codec encode input.png compressed.pde -s 0.5 -t 8.0 -q 128`

- Decode:
  - `./pde_codec decode compressed.pde output_color.png -l 10`

---

## Appendix A — Implementation mapping (functions)

- Presmoothing:
  - `gaussian_blur_2d()`, `gaussian_blur_1d()`

- Edge mask:
  - `detect_edges_marr()`

- Background anchors:
  - `add_background_samples()`
  - `detect_uniform_regions()`

- Mask compression:
  - `encode_edge_mask_jbig()`, `decode_edge_mask_jbig()`

- Signal:
  - `extract_color_signal()`, `scatter_color_signal()`

- Quantization:
  - `max_lloyd_quantization()`, `quantize_channel()`

- PDE solver:
  - `relax()`, `restrict_grid()`, `wCycle()`, `prolongate()`

- Codec:
  - `pde_encode()`, `pde_decode()`
