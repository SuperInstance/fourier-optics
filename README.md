# Fourier Optics

**A Rust library for Fourier optics computations**, including the Discrete Fourier Transform (DFT), inverse DFT, 1D Fresnel diffraction propagation, and power spectrum analysis.

## Why It Matters

Fourier optics is the mathematical foundation of optical engineering — lens design, holography, diffraction pattern analysis, telescope optics, and laser beam propagation all rely on the Fourier transform relationship between spatial fields and angular spectra. In modern applications, computational imaging (digital holography, phase retrieval, structured illumination microscopy) requires numerical implementations of these transforms. This library provides those primitives in pure Rust with no dependencies, making it suitable for embedded optical instrumentation and real-time beam profiling systems.

## How It Works

The DFT implementation follows the direct definition: for each frequency bin `k`, it sums `signal[t] · exp(-2πi·k·t/N)` across all time samples. This is a naive **O(N²)** algorithm (not FFT), which is pedagogically clear and useful for small signal sizes. The inverse DFT reverses this with a positive exponent and divides by N. Round-trip fidelity is verified in tests.

The Fresnel diffraction propagator implements the Rayleigh-Sommerfeld integral numerically: for a 1D field sampled at spacing `dx`, it computes the field at distance `z` by summing spherical wave contributions `exp(ikr)/r` from each input point to each output point, where `r = √((x_out - x_in)² + z²)` and `k = 2π/λ`. This is **O(N²)** per propagation step — accurate but computationally intensive for large apertures.

The power spectrum is simply `|F(ω)|² = Re² + Im²`, giving the energy at each frequency.

## Quick Start

```rust
use fourier_optics::{dft, idft, power_spectrum, fresnel_propagate};

fn main() {
    // Analyze a signal's frequency content
    let signal = vec![1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0];
    let spectrum = dft(&signal);

    // Power spectrum
    let power = power_spectrum(&spectrum);
    println!("DC component power: {}", power[0]);

    // Fresnel diffraction: propagate a uniform field
    let field = vec![1.0; 64];
    let propagated = fresnel_propagate(&field, 632.8e-9, 0.1, 1e-6);
    println!("Propagated field: {} samples", propagated.len());

    // Verify round-trip DFT
    let recovered = idft(&spectrum);
    for (a, b) in signal.iter().zip(&recovered) {
        assert!((a - b).abs() < 1e-10);
    }
}
```

## API

| Function | Complexity | Description |
|---|---|---|
| `dft(signal)` | **O(N²)** | Forward DFT returning `(real, imag)` pairs |
| `idft(spectrum)` | **O(N²)** | Inverse DFT recovering the original signal |
| `fresnel_propagate(field, λ, z, dx)` | **O(N²)** | 1D Fresnel diffraction propagation |
| `power_spectrum(spectrum)` | **O(N)** | Magnitude-squared of the frequency spectrum |

## Architecture Notes

Part of the SuperInstance physics simulation suite. Companion crates include `gauss-markov` (stochastic processes) and `gibbs-phenomenon` (Fourier series). See the [Architecture Guide](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT
