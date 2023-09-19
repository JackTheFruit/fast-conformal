# fast-conformal
This program computes a spherical conformal mapping for the given mesh input. The program utilizes the fast optimization technique using spherical constraints given in the paper "Folding-Free Global Conformal Mapping for Genus-0 Surfaces by Harmonic Energy Minimization" by Lai et. al (2014) and also incorporates Barzilai–Borwein step size computation for further speed.

Lai, R., Wen, Z., Yin, W., Gu, X., & Lui, L. M. (2014). Folding-free global conformal mapping for genus-0 surfaces by harmonic energy minimization. Journal of Scientific Computing, 58, 705-725.

## Setup and Usage

Run the following commands in your bash terminal

1. cd ./Code/Conformal/
2. make
3. Conformal <input_file_path> <output_file_path>

Both input and output files should be .obj or .m file formats.