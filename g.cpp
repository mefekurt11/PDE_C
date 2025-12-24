#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// Structure to represent an image grid at any level of the hierarchy
struct Grid {
    int w, h;
    std::vector<float> u;    // The solution (inpainted image)
    std::vector<float> f;    // The known pixel values (Dirichlet boundaries)
    std::vector<float> mask; // Binary mask c(x): 1 if known, 0 if unknown
    std::vector<float> r;    // Residual for multigrid cycles

    Grid(int width, int height) : w(width), h(height), 
        u(width * height, 0.0f), f(width * height, 0.0f), 
        mask(width * height, 0.0f), r(width * height, 0.0f) {}
};

// 1. Gauß-Seidel Relaxation (The basic solver) [cite: 454, 455]
void relax(Grid& g, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        for (int y = 0; y < g.h; ++y) {
            for (int x = 0; x < g.w; ++x) {
                int i = y * g.w + x;
                
                // If c(x) = 1, the value is known and "recomputed" to f [cite: 350, 461]
                if (g.mask[i] > 0.5f) {
                    g.u[i] = g.f[i];
                } else {
                    // Discrete Laplace equation: Delta u = 0 [cite: 328, 462]
                    // Uses 4-neighbourhood with reflecting Neumann boundaries 
                    float sum = 0;
                    int count = 0;
                    if (x > 0) { sum += g.u[i - 1]; count++; }
                    if (x < g.w - 1) { sum += g.u[i + 1]; count++; }
                    if (y > 0) { sum += g.u[i - g.w]; count++; }
                    if (y < g.h - 1) { sum += g.u[i + g.w]; count++; }
                    
                    if (count > 0) g.u[i] = sum / (float)count;
                }
            }
        }
    }
}

// 2. Restriction with Normalization for Sparse Data [cite: 513, 514, 517]
Grid restrict(const Grid& fine) {
    Grid coarse(fine.w / 2, fine.h / 2);
    for (int y = 0; y < coarse.h; ++y) {
        for (int x = 0; x < coarse.w; ++x) {
            int ci = y * coarse.w + x;
            int fi = (2 * y) * fine.w + (2 * x);

            // Simple area-based averaging for the mask [cite: 512]
            float m_sum = fine.mask[fi] + fine.mask[fi+1] + 
                          fine.mask[fi+fine.w] + fine.mask[fi+fine.w+1];
            coarse.mask[ci] = (m_sum > 0.1f) ? 1.0f : 0.0f; // Binary again [cite: 519]

            // Normalized convolution for the image data 
            if (m_sum > 0.0f) {
                float f_sum = fine.f[fi]*fine.mask[fi] + fine.f[fi+1]*fine.mask[fi+1] + 
                              fine.f[fi+fine.w]*fine.mask[fi+fine.w] + 
                              fine.f[fi+fine.w+1]*fine.mask[fi+fine.w+1];
                coarse.f[ci] = f_sum / m_sum;
            }
        }
    }
    return coarse;
}

// 3. Prolongation (Bilinear Interpolation) [cite: 512, 520]
void prolongateAndCorrect(const Grid& coarse, Grid& fine) {
    for (int y = 0; y < fine.h; ++y) {
        for (int x = 0; x < fine.w; ++x) {
            int fi = y * fine.w + x;
            // Only update unknown pixels [cite: 520]
            if (fine.mask[fi] < 0.5f) {
                float cx = std::min((float)x / 2.0f, (float)coarse.w - 1);
                float cy = std::min((float)y / 2.0f, (float)coarse.h - 1);
                
                // Simplified nearest-neighbor for brevity, bilinear is preferred [cite: 512]
                fine.u[fi] = coarse.u[(int)cy * coarse.w + (int)cx];
            } else {
                fine.u[fi] = fine.f[fi]; // Keep Dirichlet boundaries [cite: 520]
            }
        }
    }
}

// 4. Multigrid W-Cycle [cite: 483, 512]
void wCycle(Grid& g, int level) {
    if (g.w <= 8 || g.h <= 8 || level == 0) {
        relax(g, 50); // Solve directly on coarsest level [cite: 474]
        return;
    }

    relax(g, 2); // Presmoothing [cite: 490, 512]
    
    Grid coarse = restrict(g);
    wCycle(coarse, level - 1); // First recursive call
    wCycle(coarse, level - 1); // Second call for W-cycle [cite: 483]
    
    prolongateAndCorrect(coarse, g);
    relax(g, 2); // Postsmoothing [cite: 477, 512]
}