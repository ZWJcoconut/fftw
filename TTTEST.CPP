#include <iostream>
#include <fftw3.h>

int main() {
    // Initialize FFTW3 library
    fftw_init_threads();

    // Define input and output arrays
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 1024);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 1024);

    // Create FFTW3 plan
    fftw_plan plan = fftw_plan_dft_1d(1024, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute FFTW3 plan
    fftw_execute(plan);

    // Cleanup FFTW3 library
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup_threads();

    return 0;
}
