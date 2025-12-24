#ifndef PDE_CODEC_H
#define PDE_CODEC_H

#include <stdlib.h>
#include <math.h>
#include "jbig.h"

// Forward declarations for stb types - actual implementation in main.c
// stbi_load and stbi_write_png are used in pde_codec.c

typedef struct {
    int w, h;
    float *u;    // Reconstructed image values
    float *f;    // Known Dirichlet boundary values
    unsigned char *mask; // Binary mask c(x)
} Grid;

// ============================================================================
// A. PRESMOOTHING (Gaussian Filter) - Section 2.1
// The paper notes presmoothing is indispensable for high-quality reconstructions
// ============================================================================
void gaussian_blur_2d(int w, int h, float *image, float sigma);
void gaussian_blur_1d(float *signal, int length, float sigma);

// ============================================================================
// ENCODER FUNCTIONS
// ============================================================================

// Grid management
Grid* create_grid(int w, int h);
void free_grid(Grid *g);

// Section 2.1: Multichannel Edge Detection (Marr-Hildreth with Laplacian zero-crossings)
void detect_edges_marr(int w, int h, float *r, float *g, float *b, float threshold, unsigned char *mask);

// Background Sampling: Add Dirichlet anchor points in uniform regions
// Prevents color bleeding in large uniform backgrounds
void add_background_samples(int w, int h, const float *r, const float *g, const float *b,
                           unsigned char *mask, int sample_spacing);
void detect_uniform_regions(int w, int h, const float *r, const float *g, const float *b,
                           unsigned char *mask, float variance_threshold, int block_size);

// Section 2.2: JBIG Encoding/Decoding for binary edge mask
unsigned char* encode_edge_mask_jbig(int w, int h, const unsigned char *mask, size_t *out_len);
unsigned char* decode_edge_mask_jbig(const unsigned char *data, size_t data_len, int *out_w, int *out_h);

// Section 2.3: 1-D Signal Extraction along contours
float* extract_color_signal(int w, int h, const float *image, const unsigned char *edge_mask, 
                            int *signal_len, int **positions);
void scatter_color_signal(int w, int h, float *image, const unsigned char *edge_mask,
                          const float *signal, int signal_len, const int *positions);

// Section 2.4: Max-Lloyd Non-Uniform Quantization
void quantize_channel(float *signal, int length, int q_levels);

// ============================================================================
// DECODER FUNCTIONS (PDE-Based Inpainting via Multigrid)
// ============================================================================

// Section 5.1: Gauss-Seidel Relaxation
void relax(Grid *g, int iterations);

// Section 5.4: Normalized Restriction (coarsening)
Grid* restrict_grid(Grid *fine);

// Section 5.3: Bilinear Prolongation (interpolation)
void prolongate(Grid *coarse, Grid *fine);

// Section 5.2: Full Multigrid W-Cycle
void wCycle(Grid *g, int level);

// ============================================================================
// FULL CODEC PIPELINE
// ============================================================================

// Encode image to compressed file format
// Returns 0 on success, -1 on error
int pde_encode(const char *input_png, const char *output_file,
               float presmooth_sigma, float edge_threshold, int quant_levels);

// Decode compressed file to PNG
// Returns 0 on success, -1 on error
int pde_decode(const char *input_file, const char *output_png, int mg_levels);

#endif