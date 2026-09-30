/// 1D Discrete Fourier Transform (naive O(n²))
pub fn dft(signal: &[f64]) -> Vec<(f64, f64)> {
    let n = signal.len();
    let mut result = Vec::with_capacity(n);
    for k in 0..n {
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;
        for t in 0..n {
            let angle = -2.0 * std::f64::consts::PI * (k as f64) * (t as f64) / (n as f64);
            re += signal[t] * angle.cos();
            im += signal[t] * angle.sin();
        }
        result.push((re, im));
    }
    result
}

/// Inverse DFT
pub fn idft(spectrum: &[(f64, f64)]) -> Vec<f64> {
    let n = spectrum.len();
    let mut result = Vec::with_capacity(n);
    for t in 0..n {
        let mut val = 0.0_f64;
        for k in 0..n {
            let angle = 2.0 * std::f64::consts::PI * (k as f64) * (t as f64) / (n as f64);
            val += spectrum[k].0 * angle.cos() - spectrum[k].1 * angle.sin();
        }
        result.push(val / n as f64);
    }
    result
}

/// Fresnel diffraction integral (1D, numerical)
pub fn fresnel_propagate(field: &[f64], wavelength: f64, distance: f64, dx: f64) -> Vec<f64> {
    let n = field.len();
    let k = 2.0 * std::f64::consts::PI / wavelength;
    let mut output = vec![0.0; n];
    let z = distance;
    for i in 0..n {
        let mut sum = 0.0_f64;
        let x_out = (i as f64 - n as f64 / 2.0) * dx;
        for j in 0..n {
            let x_in = (j as f64 - n as f64 / 2.0) * dx;
            let r = ((x_out - x_in).powi(2) + z * z).sqrt();
            sum += field[j] * (k * r).cos() / r;
        }
        output[i] = sum * dx;
    }
    output
}

/// Compute power spectrum (magnitude squared)
pub fn power_spectrum(spectrum: &[(f64, f64)]) -> Vec<f64> {
    spectrum.iter().map(|(re, im)| re * re + im * im).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_dft_roundtrip() {
        let signal = vec![1.0, 2.0, 3.0, 4.0];
        let spectrum = dft(&signal);
        let recovered = idft(&spectrum);
        for (a, b) in signal.iter().zip(&recovered) {
            assert!((a - b).abs() < 1e-10);
        }
    }
}

/// FNV-1a 64 — the digest every substrate in the SuperInstance fleet agrees on.
pub const FNV_OFFSET: u64 = 0xcbf29ce484222325;
pub const FNV_PRIME: u64 = 0x100000001b3;

#[inline]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = FNV_OFFSET;
    for &b in bytes {
        h = (h ^ b as u64).wrapping_mul(FNV_PRIME);
    }
    h
}

/// True if this crate's FNV-1a still agrees with the rest of the fleet.
pub fn canary_holds() -> bool {
    fnv1a64("café Δ 日本語".as_bytes()) == 0x024a555471370b18d
}
