# fourier-optics: DFT, Inverse DFT, and Fresnel Diffraction in Rust

A zero-dependency implementation of the Discrete Fourier Transform, its inverse, Fresnel diffraction propagation, and power spectrum analysis. Designed for scalar diffraction theory simulations in one dimension.

## Why It Matters

Fourier optics is the mathematical foundation of modern imaging, holography, optical computing, and electromagnetic wave propagation analysis. The DFT connects the spatial domain (what you see) to the frequency domain (what comprises it). Fresnel diffraction predicts how a coherent wavefield evolves as it propagates through free space — critical for lens design, hologram computation, and optical interconnects in photonic computing.

## How It Works

### Discrete Fourier Transform (DFT)

The DFT transforms an N-point signal {x[t]} into N complex frequency coefficients X[k]:

```
X[k] = Σ_{t=0}^{N-1} x[t] · e^{-i·2π·k·t / N}
```

This implementation computes real and imaginary parts separately:

```
Re(X[k]) = Σ x[t] · cos(2πkt/N)
Im(X[k]) = Σ x[t] · sin(2πkt/N)
```

**Complexity**: O(N²) — the naive DFT. For production use, replace with FFT (O(N log N)).

### Inverse DFT (IDFT)

Recovers the original signal from its spectrum:

```
x[t] = (1/N) Σ_{k=0}^{N-1} [Re(X[k])·cos(2πkt/N) − Im(X[k])·sin(2πkt/N)]
```

The round-trip DFT → IDFT is exact to machine precision (tested: error < 10⁻¹⁰).

### Fresnel Diffraction

The Fresnel integral propagates a scalar field U(x') through distance z:

```
U(x) ∝ ∫ U(x') · cos(k·r) / r · dx'
```

where r = √((x−x')² + z²) and k = 2π/λ. This is the **Rayleigh–Sommerfeld** form of the scalar diffraction integral, valid for both near-field (Fresnel) and far-field (Fraunhofer) regimes.

**Complexity**: O(N²) — for each of N output points, we sum over N input points.

### Power Spectrum

```
P[k] = |X[k]|² = Re(X[k])² + Im(X[k])²
```

This is the squared magnitude of each frequency component, representing energy distribution across spatial frequencies.

### Summary Table

| Function | Time | Space | Notes |
|----------|------|-------|-------|
| `dft` | O(N²) | O(N) | Naive, no FFT |
| `idft` | O(N²) | O(N) | Exact inverse |
| `fresnel_propagate` | O(N²) | O(N) | Scalar diffraction |
| `power_spectrum` | O(N) | O(N) | Magnitude squared |

## Quick Start

```rust
use fourier_optics::{dft, idft, power_spectrum};

let signal = vec![1.0, 2.0, 3.0, 4.0];
let spectrum = dft(&signal);
let recovered = idft(&spectrum);
let power = power_spectrum(&spectrum);

// Round-trip preserves the signal
assert!((recovered[0] - 1.0).abs() < 1e-10);
```

### Fresnel Propagation

```rust
use fourier_optics::fresnel_propagate;

let aperture = vec![1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
let propagated = fresnel_propagate(&aperture, 632.8e-9, 0.01, 1e-5);
// λ = 632.8 nm (HeNe laser), z = 1 cm, dx = 10 μm
```

## API

| Function | Signature | Description |
|----------|-----------|-------------|
| `dft` | `(&[f64]) -> Vec<(f64, f64)>` | Forward DFT, returns (re, im) pairs |
| `idft` | `(&[(f64, f64)]) -> Vec<f64>` | Inverse DFT from spectrum to signal |
| `fresnel_propagate` | `(&[f64], wavelength, distance, dx) -> Vec<f64>` | Scalar Fresnel diffraction |
| `power_spectrum` | `(&[(f64, f64)]) -> Vec<f64>` | Squared magnitude per frequency |

## Architecture Notes

This crate is a **γ (gamma)** module in the γ + η = C framework — pure mathematical transforms with no side effects, no I/O, and full determinism. It provides the physical modeling layer upon which **η (eta)** orchestration (e.g., aperture optimization, beam steering) can be built. The O(N²) complexity of the naive DFT is intentionally educational; a production system would swap in `rustfft` for O(N log N) FFT computation.

## References

- Goodman, J. W. (2005). *Introduction to Fourier Optics* (3rd ed.). Roberts & Company.
- Bracewell, R. N. (2000). *The Fourier Transform and Its Applications* (3rd ed.). McGraw-Hill.
- Born, M. & Wolf, E. (2019). *Principles of Optics* (7th ed.). Cambridge University Press.

## License

MIT
