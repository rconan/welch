//! # Welch periodogram
//!
//! Welch method:
//!  - split signal in overlapping segment
//!  - DFT of segments
//!  - average segments DFT
//!
//! Properties:
//!  - \# of segments: k
//!  - overlapping factor: a
//!  - segment size: l
//!
//! signal length: n = kl - (k-1)la = l(k(1-a)+a), so l = n/(k(1-a)+a)

use num_complex::Complex;
use rustfft::FftPlanner;

pub struct Welch<'a> {
    pub n_segment: usize,
    pub overlap: f64,
    pub segment_size: usize,
    pub signal: &'a Vec<f64>,
}

impl<'a> Welch<'a> {
    pub fn new(n_segment: usize, overlap: f64, signal: &'a Vec<f64>) -> Self {
        let l =
            (signal.len() as f64 / (n_segment as f64 * (1. - overlap) + overlap)).trunc() as usize;
        Self {
            n_segment,
            overlap,
            segment_size: l,
            signal,
        }
    }
    pub fn segmenting(&self) -> Vec<Complex<f64>> {
        let l = self.segment_size;
        let a = self.overlap;
        let nel = l - (l as f64 * a).round() as usize;
        self.signal
            .windows(self.segment_size)
            .step_by(nel)
            .map(|segment| {
                segment
                    .iter()
                    .map(|x| Complex::new(*x, 0f64))
                    .collect::<Vec<Complex<f64>>>()
            })
            .flatten()
            .collect()
    }
    pub fn dft(&self) -> Vec<Complex<f64>> {
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(self.segment_size);
        let mut buffer = self.segmenting();
        fft.process(&mut buffer);
        buffer
    }
    pub fn periogram(&self) -> Vec<f64> {
        let buffer = self.dft();
        let n = self.segment_size / 2;
        /*
                let psd: Vec<_> = buffer
                    .chunks(self.segment_size)
                    .map(|x| x.iter().take(n).map(|x| x.norm_sqr()).collect::<Vec<f64>>())
                    .collect();
        */
        buffer
            .chunks(self.segment_size)
            .fold(vec![0f64; n], |a, x| {
                a.iter() // 0 0 0,  a b c
                    .zip(x.iter()) // a b c, d e f
                    .map(|(a, x)| a + x.norm_sqr())
                    .collect::<Vec<f64>>() // 0+a 0+b 0+c, a+d b+e c+f
            })
        /*
        a    : 0 0 0  0+a 0+b 0+c  a+d b+e c+f
        x[0] : a b c
        x[1] : d e f
         */
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
