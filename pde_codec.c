#include "pde_codec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Include stb headers (implementation is in main.c)
#include "stb_image.h"
#include "stb_image_write.h"

// ============================================================================
// A. PRESMOOTHING (Gaussian Filter) - Section 2.1
// The paper notes presmoothing is indispensable for high-quality reconstructions
// ============================================================================

static int clampi(int v, int lo, int hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

// Generate 1D Gaussian kernel with given sigma
static float* make_gaussian_kernel(float sigma, int *out_radius) {
    int radius = (int)ceilf(3.0f * sigma);
    if (radius < 1) radius = 1;
    int size = 2 * radius + 1;
    
    float *kernel = (float*)malloc(size * sizeof(float));
    float sum = 0.0f;
    float two_sigma2 = 2.0f * sigma * sigma;
    
    for (int i = -radius; i <= radius; ++i) {
        float v = expf(-(float)(i * i) / two_sigma2);
        kernel[i + radius] = v;
        sum += v;
    }
    // Normalize
    for (int i = 0; i < size; ++i) kernel[i] /= sum;
    
    *out_radius = radius;
    return kernel;
}

// Separable 2D Gaussian blur (in-place)
void gaussian_blur_2d(int w, int h, float *image, float sigma) {
    if (sigma <= 0.0f) return;
    
    int radius;
    float *kernel = make_gaussian_kernel(sigma, &radius);
    
    float *tmp = (float*)malloc(w * h * sizeof(float));
    
    // Horizontal pass
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            float acc = 0.0f;
            for (int k = -radius; k <= radius; ++k) {
                int xx = clampi(x + k, 0, w - 1);
                acc += image[y * w + xx] * kernel[k + radius];
            }
            tmp[y * w + x] = acc;
        }
    }
    
    // Vertical pass
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            float acc = 0.0f;
            for (int k = -radius; k <= radius; ++k) {
                int yy = clampi(y + k, 0, h - 1);
                acc += tmp[yy * w + x] * kernel[k + radius];
            }
            image[y * w + x] = acc;
        }
    }
    
    free(tmp);
    free(kernel);
}

// 1D Gaussian blur for color signal (in-place)
void gaussian_blur_1d(float *signal, int length, float sigma) {
    if (sigma <= 0.0f || length <= 1) return;
    
    int radius;
    float *kernel = make_gaussian_kernel(sigma, &radius);
    
    float *tmp = (float*)malloc(length * sizeof(float));
    
    for (int i = 0; i < length; ++i) {
        float acc = 0.0f;
        for (int k = -radius; k <= radius; ++k) {
            int ii = clampi(i + k, 0, length - 1);
            acc += signal[ii] * kernel[k + radius];
        }
        tmp[i] = acc;
    }
    
    memcpy(signal, tmp, length * sizeof(float));
    free(tmp);
    free(kernel);
}

// ============================================================================
// GRID MANAGEMENT
// ============================================================================

Grid* create_grid(int w, int h) {
    Grid *g = (Grid*)malloc(sizeof(Grid));
    g->w = w; g->h = h;
    g->u = (float*)calloc(w * h, sizeof(float));
    g->f = (float*)calloc(w * h, sizeof(float));
    g->mask = (unsigned char*)calloc(w * h, sizeof(unsigned char));
    return g;
}

void free_grid(Grid *g) {
    if (!g) return;
    free(g->u); free(g->f); free(g->mask); free(g);
}

// ============================================================================
// SECTION 2.1: MULTICHANNEL EDGE DETECTION (Marr-Hildreth)
// Uses Laplacian zero-crossings with optional presmoothing
// ============================================================================

void detect_edges_marr(int w, int h, float *r, float *g, float *b, float threshold, unsigned char *mask) {
    float *lap = (float*)calloc(w * h, sizeof(float));
    
    // Discrete Laplacian calculation for 3 channels combined (Equation 1)
    // L(u) = u_xx + u_yy approximated by 5-point stencil
    for (int y = 1; y < h - 1; ++y) {
        for (int x = 1; x < w - 1; ++x) {
            int i = y * w + x;
            // Laplacian for each channel
            float lr = r[i+1] + r[i-1] + r[i+w] + r[i-w] - 4.0f * r[i];
            float lg = g[i+1] + g[i-1] + g[i+w] + g[i-w] - 4.0f * g[i];
            float lb = b[i+1] + b[i-1] + b[i+w] + b[i-w] - 4.0f * b[i];
            // Combine: sqrt(lr^2 + lg^2 + lb^2) for magnitude, but paper uses sum
            lap[i] = sqrtf(lr*lr + lg*lg + lb*lb);
        }
    }

    // Zero-crossing detection with threshold
    for (int y = 1; y < h - 1; ++y) {
        for (int x = 1; x < w - 1; ++x) {
            int i = y * w + x;
            
            // Edge if gradient magnitude exceeds threshold
            if (lap[i] > threshold) {
                mask[i] = 1;
            }
        }
    }
    free(lap);
}

// ============================================================================
// BACKGROUND SAMPLING: Add Dirichlet anchor points in uniform regions
// The paper treats background as "missing data" filled by Laplace equation.
// For large uniform regions, we need anchor points to prevent color bleeding.
// ============================================================================

// Add regular grid samples in regions far from any edge
void add_background_samples(int w, int h, const float *r, const float *g, const float *b,
                           unsigned char *mask, int sample_spacing) {
    if (sample_spacing <= 0) return;
    
    // For each potential sample point
    for (int y = sample_spacing; y < h - sample_spacing; y += sample_spacing) {
        for (int x = sample_spacing; x < w - sample_spacing; x += sample_spacing) {
            int i = y * w + x;
            
            // Skip if already an edge
            if (mask[i]) continue;
            
            // Check if there's an edge within sample_spacing/2 radius
            int has_nearby_edge = 0;
            int check_radius = sample_spacing / 2;
            
            for (int dy = -check_radius; dy <= check_radius && !has_nearby_edge; dy++) {
                for (int dx = -check_radius; dx <= check_radius && !has_nearby_edge; dx++) {
                    int ny = y + dy;
                    int nx = x + dx;
                    if (ny >= 0 && ny < h && nx >= 0 && nx < w) {
                        if (mask[ny * w + nx]) {
                            has_nearby_edge = 1;
                        }
                    }
                }
            }
            
            // If no nearby edge, add this as an anchor point
            if (!has_nearby_edge) {
                mask[i] = 1;
            }
        }
    }
}

// Detect uniform color regions and add sample points
// This ensures large flat backgrounds are properly anchored
void detect_uniform_regions(int w, int h, const float *r, const float *g, const float *b,
                           unsigned char *mask, float variance_threshold, int block_size) {
    // Process image in blocks
    for (int by = 0; by < h; by += block_size) {
        for (int bx = 0; bx < w; bx += block_size) {
            // Calculate color variance in this block
            float sum_r = 0, sum_g = 0, sum_b = 0;
            float sum_r2 = 0, sum_g2 = 0, sum_b2 = 0;
            int count = 0;
            int has_edge = 0;
            
            for (int y = by; y < by + block_size && y < h; y++) {
                for (int x = bx; x < bx + block_size && x < w; x++) {
                    int i = y * w + x;
                    if (mask[i]) has_edge = 1;
                    sum_r += r[i]; sum_g += g[i]; sum_b += b[i];
                    sum_r2 += r[i]*r[i]; sum_g2 += g[i]*g[i]; sum_b2 += b[i]*b[i];
                    count++;
                }
            }
            
            // Skip blocks that already have edges
            if (has_edge || count == 0) continue;
            
            float mean_r = sum_r / count;
            float mean_g = sum_g / count;
            float mean_b = sum_b / count;
            
            float var_r = (sum_r2 / count) - (mean_r * mean_r);
            float var_g = (sum_g2 / count) - (mean_g * mean_g);
            float var_b = (sum_b2 / count) - (mean_b * mean_b);
            
            float total_var = var_r + var_g + var_b;
            
            // If block is uniform (low variance), add center point as Dirichlet anchor
            if (total_var < variance_threshold) {
                int cx = bx + block_size / 2;
                int cy = by + block_size / 2;
                if (cx < w && cy < h) {
                    mask[cy * w + cx] = 1;
                }
            }
        }
    }
}

// ============================================================================
// SECTION 2.2: JBIG ENCODING/DECODING FOR BINARY EDGE MASK
// Paper recommends JBIG for lossless compression of binary edge images
// ============================================================================

// Callback for JBIG encoder output
static void jbig_output_callback(unsigned char *data, size_t len, void *user) {
    // Append to dynamic buffer
    typedef struct { unsigned char *buf; size_t len; size_t cap; } Buffer;
    Buffer *b = (Buffer*)user;
    while (b->len + len > b->cap) {
        b->cap = b->cap ? b->cap * 2 : 4096;
        b->buf = (unsigned char*)realloc(b->buf, b->cap);
    }
    memcpy(b->buf + b->len, data, len);
    b->len += len;
}

unsigned char* encode_edge_mask_jbig(int w, int h, const unsigned char *mask, size_t *out_len) {
    // Convert mask to packed bitmap (1 bit per pixel, row-padded to bytes)
    int bytes_per_row = (w + 7) / 8;
    unsigned char *bitmap = (unsigned char*)calloc(bytes_per_row * h, 1);
    
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            if (mask[y * w + x]) {
                int byte_idx = y * bytes_per_row + x / 8;
                int bit_idx = 7 - (x % 8);
                bitmap[byte_idx] |= (1 << bit_idx);
            }
        }
    }
    
    // Setup JBIG encoder
    struct jbg_enc_state state;
    unsigned char *bptr = bitmap;
    
    typedef struct { unsigned char *buf; size_t len; size_t cap; } Buffer;
    Buffer buffer = {NULL, 0, 0};
    
    jbg_enc_init(&state, w, h, 1, &bptr, jbig_output_callback, &buffer);
    
    // Configure: single resolution layer, typical settings for edge images
    jbg_enc_options(&state, JBG_ILEAVE | JBG_SMID, 0, -1, -1, -1);
    
    // Encode
    jbg_enc_out(&state);
    jbg_enc_free(&state);
    free(bitmap);
    
    *out_len = buffer.len;
    return buffer.buf;
}

unsigned char* decode_edge_mask_jbig(const unsigned char *data, size_t data_len, int *out_w, int *out_h) {
    struct jbg_dec_state state;
    jbg_dec_init(&state);
    
    // Cast away const - jbig API doesn't modify the data but has incorrect signature
    int result = jbg_dec_in(&state, (unsigned char*)data, data_len, NULL);
    if (result != JBG_EOK && result != JBG_EOK_INTR) {
        jbg_dec_free(&state);
        return NULL;
    }
    
    int w = jbg_dec_getwidth(&state);
    int h = jbg_dec_getheight(&state);
    unsigned char *bitmap = jbg_dec_getimage(&state, 0);
    int bytes_per_row = (w + 7) / 8;
    
    // Convert packed bitmap back to byte-per-pixel mask
    unsigned char *mask = (unsigned char*)calloc(w * h, 1);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int byte_idx = y * bytes_per_row + x / 8;
            int bit_idx = 7 - (x % 8);
            mask[y * w + x] = (bitmap[byte_idx] >> bit_idx) & 1;
        }
    }
    
    *out_w = w;
    *out_h = h;
    jbg_dec_free(&state);
    return mask;
}

// ============================================================================
// SECTION 2.3: 1-D SIGNAL EXTRACTION ALONG CONTOURS
// Extract color values from edge pixels in raster order for entropy coding
// ============================================================================

float* extract_color_signal(int w, int h, const float *image, const unsigned char *edge_mask,
                            int *signal_len, int **positions) {
    // Count edge + border pixels
    int count = 0;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int i = y * w + x;
            int is_border = (x == 0 || x == w-1 || y == 0 || y == h-1);
            if (edge_mask[i] || is_border) count++;
        }
    }
    
    float *signal = (float*)malloc(count * sizeof(float));
    int *pos = (int*)malloc(count * sizeof(int));
    
    int idx = 0;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int i = y * w + x;
            int is_border = (x == 0 || x == w-1 || y == 0 || y == h-1);
            if (edge_mask[i] || is_border) {
                signal[idx] = image[i];
                pos[idx] = i;
                idx++;
            }
        }
    }
    
    *signal_len = count;
    if (positions) *positions = pos;
    else free(pos);
    
    return signal;
}

void scatter_color_signal(int w, int h, float *image, const unsigned char *edge_mask,
                          const float *signal, int signal_len, const int *positions) {
    for (int i = 0; i < signal_len; ++i) {
        image[positions[i]] = signal[i];
    }
}

// ============================================================================
// SECTION 2.4: MAX-LLOYD NON-UNIFORM QUANTIZATION
// ============================================================================

void max_lloyd_quantization(float *signal, int length, int q_levels, float *recon_points) {
    // Initial uniform distribution
    for (int i = 0; i < q_levels; ++i) recon_points[i] = (255.0f / (q_levels - 1)) * i;

    float boundaries[256];
    for (int iter = 0; iter < 100; ++iter) {
        for (int i = 0; i < q_levels - 1; ++i) 
            boundaries[i] = (recon_points[i] + recon_points[i+1]) / 2.0f;

        float sums[256] = {0};
        int counts[256] = {0};

        for (int i = 0; i < length; ++i) {
            int idx = 0;
            while (idx < q_levels - 1 && signal[i] > boundaries[idx]) idx++;
            sums[idx] += signal[i];
            counts[idx]++;
        }

        float max_shift = 0;
        for (int i = 0; i < q_levels; ++i) {
            if (counts[i] > 0) {
                float new_pt = sums[i] / counts[i];
                max_shift = fmaxf(max_shift, fabsf(new_pt - recon_points[i]));
                recon_points[i] = new_pt;
            }
        }
        if (max_shift < 0.01f) break;
    }
}

void quantize_channel(float *signal, int length, int q_levels) {
    if (q_levels <= 1) return;
    float recon_pts[256];
    max_lloyd_quantization(signal, length, q_levels, recon_pts);
    
    float boundaries[256];
    for (int i = 0; i < q_levels - 1; ++i) boundaries[i] = (recon_pts[i] + recon_pts[i+1]) / 2.0f;

    for (int i = 0; i < length; ++i) {
        int idx = 0;
        while (idx < q_levels - 1 && signal[i] > boundaries[idx]) idx++;
        signal[i] = recon_pts[idx];
    }
}

// --- SECTION 5: MULTIGRID SOLVER (HOMOGENEOUS DIFFUSION) ---

// Section 5.4: Normalized Restriction (Equation 22)
Grid* restrict_grid(Grid *fine) {
    Grid *coarse = create_grid(fine->w / 2, fine->h / 2);
    for (int y = 0; y < coarse->h; ++y) {
        for (int x = 0; x < coarse->w; ++x) {
            int ci = y * coarse->w + x;
            double f_sum = 0, m_sum = 0, u_sum = 0;
            for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                    int fi = (2 * y + dy) * fine->w + (2 * x + dx);
                    if (fine->mask[fi]) {
                        f_sum += fine->f[fi];
                        m_sum += 1.0;
                    }
                    u_sum += fine->u[fi];
                }
            }
            if (m_sum > 0) {
                coarse->mask[ci] = 1;
                coarse->f[ci] = (float)(f_sum / m_sum);
                coarse->u[ci] = coarse->f[ci]; 
            } else {
                coarse->mask[ci] = 0;
                coarse->u[ci] = (float)(u_sum / 4.0);
            }
        }
    }
    return coarse;
}

// Section 5.1: High-Precision Relaxation
void relax(Grid *g, int iterations) {
    for (int it = 0; it < iterations; ++it) {
        for (int y = 0; y < g->h; ++y) {
            for (int x = 0; x < g->w; ++x) {
                int i = y * g->w + x;
                if (g->mask[i]) {
                    g->u[i] = g->f[i];
                } else {
                    float sum = 0; int cnt = 0;
                    if (x > 0) { sum += g->u[i-1]; cnt++; }
                    if (x < g->w-1) { sum += g->u[i+1]; cnt++; }
                    if (y > 0) { sum += g->u[i-g->w]; cnt++; }
                    if (y < g->h-1) { sum += g->u[i+g->w]; cnt++; }
                    if (cnt > 0) g->u[i] = sum / (float)cnt;
                }
            }
        }
    }
}



void wCycle(Grid *g, int level) {
    if (level == 0 || g->w < 8 || g->h < 8) {
        relax(g, 500); // Heavy relaxation at coarsest level
        return;
    }
    relax(g, 3); // Pre-smoothing
    Grid *coarse = restrict_grid(g);
    wCycle(coarse, level - 1);
    wCycle(coarse, level - 1); 
    for (int i = 0; i < g->w * g->h; i++) {
        if (!g->mask[i]) g->u[i] = coarse->u[((i/g->w)/2)*coarse->w + (i%g->w)/2];
    }
    relax(g, 4);
    free_grid(coarse);
}

// ============================================================================
// SECTION 5.2: BILINEAR PROLONGATION
// Interpolates coarse grid values to fine grid
// ============================================================================

void prolongate(Grid *coarse, Grid *fine) {
    // Bilinear interpolation from coarse to fine grid
    for (int y = 0; y < fine->h; ++y) {
        for (int x = 0; x < fine->w; ++x) {
            int fi = y * fine->w + x;
            if (fine->mask[fi]) continue; // Preserve known data points
            
            // Map fine coords to coarse coords
            float cx = (float)x / 2.0f;
            float cy = (float)y / 2.0f;
            
            int x0 = (int)floorf(cx);
            int y0 = (int)floorf(cy);
            int x1 = x0 + 1;
            int y1 = y0 + 1;
            
            // Clamp to coarse grid bounds
            if (x0 < 0) x0 = 0;
            if (y0 < 0) y0 = 0;
            if (x1 >= coarse->w) x1 = coarse->w - 1;
            if (y1 >= coarse->h) y1 = coarse->h - 1;
            
            float dx = cx - x0;
            float dy = cy - y0;
            
            // Bilinear interpolation
            float v00 = coarse->u[y0 * coarse->w + x0];
            float v10 = coarse->u[y0 * coarse->w + x1];
            float v01 = coarse->u[y1 * coarse->w + x0];
            float v11 = coarse->u[y1 * coarse->w + x1];
            
            fine->u[fi] = (1-dx)*(1-dy)*v00 + dx*(1-dy)*v10 + (1-dx)*dy*v01 + dx*dy*v11;
        }
    }
}

// ============================================================================
// FULL CODEC PIPELINE
// ============================================================================

int pde_encode(const char *input_png, const char *output_file,
               float presmooth_sigma, float edge_threshold, int quant_levels) {
    // Load image
    int w, h, channels;
    unsigned char *img_data = stbi_load(input_png, &w, &h, &channels, 3);
    if (!img_data) {
        fprintf(stderr, "Error: Cannot load %s\n", input_png);
        return -1;
    }
    
    printf("Encoding %s (%dx%d)\n", input_png, w, h);
    
    // Separate channels to float (original values)
    float *R = (float*)malloc(w * h * sizeof(float));
    float *G = (float*)malloc(w * h * sizeof(float));
    float *B = (float*)malloc(w * h * sizeof(float));
    
    for (int i = 0; i < w * h; ++i) {
        R[i] = (float)img_data[i * 3 + 0];
        G[i] = (float)img_data[i * 3 + 1];
        B[i] = (float)img_data[i * 3 + 2];
    }
    stbi_image_free(img_data);
    
    // A. Presmoothing for edge detection only
    // Create smoothed copies for edge detection
    float *R_smooth = (float*)malloc(w * h * sizeof(float));
    float *G_smooth = (float*)malloc(w * h * sizeof(float));
    float *B_smooth = (float*)malloc(w * h * sizeof(float));
    memcpy(R_smooth, R, w * h * sizeof(float));
    memcpy(G_smooth, G, w * h * sizeof(float));
    memcpy(B_smooth, B, w * h * sizeof(float));
    
    if (presmooth_sigma > 0) {
        printf("  Presmoothing with sigma=%.2f\n", presmooth_sigma);
        gaussian_blur_2d(w, h, R_smooth, presmooth_sigma);
        gaussian_blur_2d(w, h, G_smooth, presmooth_sigma);
        gaussian_blur_2d(w, h, B_smooth, presmooth_sigma);
    }
    
    // B. Edge detection (Marr-Hildreth) on smoothed image
    unsigned char *edge_mask = (unsigned char*)calloc(w * h, 1);
    detect_edges_marr(w, h, R_smooth, G_smooth, B_smooth, edge_threshold, edge_mask);
    
    free(R_smooth);
    free(G_smooth);
    free(B_smooth);
    
    // Count edge pixels
    int edge_count = 0;
    for (int i = 0; i < w * h; ++i) if (edge_mask[i]) edge_count++;
    printf("  Detected %d edge pixels (%.2f%%)\n", edge_count, 100.0f * edge_count / (w * h));
    
    // B2. ADD BACKGROUND SAMPLES (Dirichlet anchor points)
    // This is crucial for uniform regions like backgrounds.
    // The Laplace equation Δu=0 fills unknown regions, but needs anchor points.
    int bg_sample_spacing = 16;  // Sample spacing for background regions (smaller = more samples)
    add_background_samples(w, h, R, G, B, edge_mask, bg_sample_spacing);
    
    // Also detect large uniform regions and add samples
    float variance_threshold = 100.0f;  // Low variance = uniform region
    int block_size = 24;  // Smaller blocks for finer sampling
    detect_uniform_regions(w, h, R, G, B, edge_mask, variance_threshold, block_size);
    
    // Recount total sample points
    int total_samples = 0;
    for (int i = 0; i < w * h; ++i) if (edge_mask[i]) total_samples++;
    printf("  Total Dirichlet points (edges + background): %d (%.2f%%)\n", 
           total_samples, 100.0f * total_samples / (w * h));
    
    // C. Encode edge mask with JBIG
    size_t jbig_len;
    unsigned char *jbig_data = encode_edge_mask_jbig(w, h, edge_mask, &jbig_len);
    printf("  JBIG edge mask: %zu bytes\n", jbig_len);
    
    // D. Extract 1-D color signals from ORIGINAL (not smoothed) image
    int sig_len;
    int *positions;
    float *sig_R = extract_color_signal(w, h, R, edge_mask, &sig_len, &positions);
    float *sig_G = extract_color_signal(w, h, G, edge_mask, &sig_len, NULL);
    float *sig_B = extract_color_signal(w, h, B, edge_mask, &sig_len, NULL);
    printf("  Signal length: %d samples\n", sig_len);
    
    // E. Optional smoothing of 1-D signals for better compression
    // Paper suggests this helps with quantization
    if (presmooth_sigma > 0) {
        gaussian_blur_1d(sig_R, sig_len, presmooth_sigma * 0.5f);
        gaussian_blur_1d(sig_G, sig_len, presmooth_sigma * 0.5f);
        gaussian_blur_1d(sig_B, sig_len, presmooth_sigma * 0.5f);
    }
    
    // F. Quantize signals with Max-Lloyd
    if (quant_levels > 0 && quant_levels < 256) {
        printf("  Quantizing to %d levels\n", quant_levels);
        quantize_channel(sig_R, sig_len, quant_levels);
        quantize_channel(sig_G, sig_len, quant_levels);
        quantize_channel(sig_B, sig_len, quant_levels);
    }
    
    // G. Write compressed file
    FILE *fp = fopen(output_file, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create %s\n", output_file);
        return -1;
    }
    
    // Header: magic, dimensions, parameters
    const char magic[] = "PDE1";
    fwrite(magic, 1, 4, fp);
    fwrite(&w, sizeof(int), 1, fp);
    fwrite(&h, sizeof(int), 1, fp);
    fwrite(&presmooth_sigma, sizeof(float), 1, fp);
    fwrite(&edge_threshold, sizeof(float), 1, fp);
    fwrite(&quant_levels, sizeof(int), 1, fp);
    
    // JBIG data
    fwrite(&jbig_len, sizeof(size_t), 1, fp);
    fwrite(jbig_data, 1, jbig_len, fp);
    
    // Signal data (as 8-bit after quantization)
    fwrite(&sig_len, sizeof(int), 1, fp);
    for (int i = 0; i < sig_len; ++i) {
        unsigned char r8 = (unsigned char)clampi((int)roundf(sig_R[i]), 0, 255);
        unsigned char g8 = (unsigned char)clampi((int)roundf(sig_G[i]), 0, 255);
        unsigned char b8 = (unsigned char)clampi((int)roundf(sig_B[i]), 0, 255);
        fwrite(&r8, sizeof(unsigned char), 1, fp);
        fwrite(&g8, sizeof(unsigned char), 1, fp);
        fwrite(&b8, sizeof(unsigned char), 1, fp);
    }
    
    long file_size = ftell(fp);
    fclose(fp);
    
    printf("  Output file: %ld bytes (%.2f bpp)\n", file_size, 8.0f * file_size / (w * h));
    
    // Cleanup
    free(R); free(G); free(B);
    free(edge_mask);
    free(jbig_data);
    free(sig_R); free(sig_G); free(sig_B);
    free(positions);
    
    return 0;
}

int pde_decode(const char *input_file, const char *output_png, int mg_levels) {
    FILE *fp = fopen(input_file, "rb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", input_file);
        return -1;
    }
    
    // Read header
    char magic[4];
    fread(magic, 1, 4, fp);
    if (memcmp(magic, "PDE1", 4) != 0) {
        fprintf(stderr, "Error: Invalid file format\n");
        fclose(fp);
        return -1;
    }
    
    int w, h;
    float presmooth_sigma, edge_threshold;
    int quant_levels;
    fread(&w, sizeof(int), 1, fp);
    fread(&h, sizeof(int), 1, fp);
    fread(&presmooth_sigma, sizeof(float), 1, fp);
    fread(&edge_threshold, sizeof(float), 1, fp);
    fread(&quant_levels, sizeof(int), 1, fp);
    
    printf("Decoding to %s (%dx%d)\n", output_png, w, h);
    
    // Read JBIG data
    size_t jbig_len;
    fread(&jbig_len, sizeof(size_t), 1, fp);
    unsigned char *jbig_data = (unsigned char*)malloc(jbig_len);
    fread(jbig_data, 1, jbig_len, fp);
    
    // Decode edge mask
    int dw, dh;
    unsigned char *edge_mask = decode_edge_mask_jbig(jbig_data, jbig_len, &dw, &dh);
    free(jbig_data);
    
    if (!edge_mask || dw != w || dh != h) {
        fprintf(stderr, "Error: JBIG decode failed\n");
        fclose(fp);
        return -1;
    }
    
    // Read signal data
    int sig_len;
    fread(&sig_len, sizeof(int), 1, fp);
    
    float *sig_R = (float*)malloc(sig_len * sizeof(float));
    float *sig_G = (float*)malloc(sig_len * sizeof(float));
    float *sig_B = (float*)malloc(sig_len * sizeof(float));
    
    for (int i = 0; i < sig_len; ++i) {
        unsigned char r8, g8, b8;
        fread(&r8, sizeof(unsigned char), 1, fp);
        fread(&g8, sizeof(unsigned char), 1, fp);
        fread(&b8, sizeof(unsigned char), 1, fp);
        sig_R[i] = (float)r8;
        sig_G[i] = (float)g8;
        sig_B[i] = (float)b8;
    }
    fclose(fp);
    
    // Reconstruct positions from edge mask
    int *positions = (int*)malloc(sig_len * sizeof(int));
    int idx = 0;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int i = y * w + x;
            int is_border = (x == 0 || x == w-1 || y == 0 || y == h-1);
            if (edge_mask[i] || is_border) {
                positions[idx++] = i;
            }
        }
    }
    
    // Create grids for each channel
    Grid *grid_R = create_grid(w, h);
    Grid *grid_G = create_grid(w, h);
    Grid *grid_B = create_grid(w, h);
    
    // Set known values (edges + borders)
    for (int i = 0; i < sig_len; ++i) {
        int pos = positions[i];
        grid_R->f[pos] = sig_R[i];
        grid_G->f[pos] = sig_G[i];
        grid_B->f[pos] = sig_B[i];
        grid_R->mask[pos] = 1;
        grid_G->mask[pos] = 1;
        grid_B->mask[pos] = 1;
        grid_R->u[pos] = sig_R[i];
        grid_G->u[pos] = sig_G[i];
        grid_B->u[pos] = sig_B[i];
    }
    
    // Run multigrid W-cycle solver (multiple cycles for convergence)
    printf("  Running W-cycle solver (%d levels)...\n", mg_levels);
    for (int cycle = 0; cycle < 3; ++cycle) {
        wCycle(grid_R, mg_levels);
        wCycle(grid_G, mg_levels);
        wCycle(grid_B, mg_levels);
    }
    // Final relaxation pass
    relax(grid_R, 50);
    relax(grid_G, 50);
    relax(grid_B, 50);
    printf("  Done.\n");
    
    // Compose output image
    unsigned char *out_img = (unsigned char*)malloc(w * h * 3);
    for (int i = 0; i < w * h; ++i) {
        out_img[i * 3 + 0] = (unsigned char)clampi((int)roundf(grid_R->u[i]), 0, 255);
        out_img[i * 3 + 1] = (unsigned char)clampi((int)roundf(grid_G->u[i]), 0, 255);
        out_img[i * 3 + 2] = (unsigned char)clampi((int)roundf(grid_B->u[i]), 0, 255);
    }
    
    stbi_write_png(output_png, w, h, 3, out_img, w * 3);
    printf("  Saved to %s\n", output_png);
    
    // Cleanup
    free_grid(grid_R);
    free_grid(grid_G);
    free_grid(grid_B);
    free(sig_R); free(sig_G); free(sig_B);
    free(positions);
    free(edge_mask);
    free(out_img);
    
    return 0;
}