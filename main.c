#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "pde_codec.h"

void print_usage(const char *prog) {
    printf("PDE Edge-Based Image Compression\n");
    printf("Based on Mainberger et al. paper\n\n");
    printf("Usage:\n");
    printf("  %s encode <input.png> <output.pde> [options]\n", prog);
    printf("  %s decode <input.pde> <output.png> [options]\n", prog);
    printf("  %s demo                             - Run demo on input.png\n\n", prog);
    printf("Encode options:\n");
    printf("  -s <sigma>      Presmoothing sigma (default: 0.8)\n");
    printf("  -t <threshold>  Edge detection threshold (default: 15.0)\n");
    printf("  -q <levels>     Quantization levels (default: 64)\n\n");
    printf("Decode options:\n");
    printf("  -l <levels>     Multigrid levels (default: 8)\n");
}

int main(int argc, char *argv[]) {
    // Default parameters (tuned according to paper)
    float presmooth_sigma = 0.5f;   // Section A: presmoothing (lower for less blur)
    float edge_threshold = 8.0f;    // Section 2.1: edge detection (balanced)
    int quant_levels = 128;         // Section 2.4: quantization (128 for good quality)
    int mg_levels = 10;             // Section 5: multigrid levels
    
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    if (strcmp(argv[1], "encode") == 0) {
        if (argc < 4) {
            printf("Error: encode requires input and output files\n");
            print_usage(argv[0]);
            return 1;
        }
        
        const char *input_file = argv[2];
        const char *output_file = argv[3];
        
        // Parse options
        for (int i = 4; i < argc; i += 2) {
            if (i + 1 >= argc) break;
            if (strcmp(argv[i], "-s") == 0) presmooth_sigma = atof(argv[i+1]);
            else if (strcmp(argv[i], "-t") == 0) edge_threshold = atof(argv[i+1]);
            else if (strcmp(argv[i], "-q") == 0) quant_levels = atoi(argv[i+1]);
        }
        
        return pde_encode(input_file, output_file, presmooth_sigma, edge_threshold, quant_levels);
        
    } else if (strcmp(argv[1], "decode") == 0) {
        if (argc < 4) {
            printf("Error: decode requires input and output files\n");
            print_usage(argv[0]);
            return 1;
        }
        
        const char *input_file = argv[2];
        const char *output_file = argv[3];
        
        // Parse options
        for (int i = 4; i < argc; i += 2) {
            if (i + 1 >= argc) break;
            if (strcmp(argv[i], "-l") == 0) mg_levels = atoi(argv[i+1]);
        }
        
        return pde_decode(input_file, output_file, mg_levels);
        
    } else if (strcmp(argv[1], "demo") == 0) {
        // Demo mode: encode and decode input.png
        printf("=== PDE Compression Demo ===\n\n");
        
        // Load original for comparison
        int w, h, ch;
        unsigned char *original = stbi_load("input.png", &w, &h, &ch, 3);
        if (!original) {
            printf("Error: Cannot load input.png\n");
            return -1;
        }
        
        // Encode
        printf("[1] Encoding...\n");
        if (pde_encode("input.png", "compressed.pde", presmooth_sigma, edge_threshold, quant_levels) != 0) {
            stbi_image_free(original);
            return -1;
        }
        
        // Decode
        printf("\n[2] Decoding...\n");
        if (pde_decode("compressed.pde", "output_color.png", mg_levels) != 0) {
            stbi_image_free(original);
            return -1;
        }
        
        // Load reconstructed and compute PSNR
        unsigned char *reconstructed = stbi_load("output_color.png", &w, &h, &ch, 3);
        if (reconstructed) {
            double mse = 0;
            for (int i = 0; i < w * h * 3; ++i) {
                double diff = (double)original[i] - (double)reconstructed[i];
                mse += diff * diff;
            }
            mse /= (3.0 * w * h);
            double psnr = 10.0 * log10((255.0 * 255.0) / mse);
            
            printf("\n=== RESULTS ===\n");
            printf("PSNR: %.2f dB\n", psnr);
            
            stbi_image_free(reconstructed);
        }
        
        stbi_image_free(original);
        printf("\nDone! Output saved to output_color.png\n");
        return 0;
        
    } else {
        printf("Unknown command: %s\n", argv[1]);
        print_usage(argv[0]);
        return 1;
    }
}