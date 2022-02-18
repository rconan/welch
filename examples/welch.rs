use std::fmt::Display;
use welch::{Welch, Window};

#[derive(Debug)]
pub struct One {
    weights: Vec<f64>,
}
impl Window for One {
    fn new(n: usize) -> Self {
        Self {
            weights: vec![1f64; n],
        }
    }
    fn weights(&self) -> &[f64] {
        self.weights.as_slice()
    }
}
impl One {
    pub fn next(&self) {}
}
impl Display for One {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "I am the One")
    }
}

fn main() {
    // sin(2\pi i f0 / sampling_frequency) i: [0 , sampling_frequency/f0]
    let sampling_frequency = 5_000f64; // 1_000
    let f0 = 150f64;
    let n_step = 10 * (sampling_frequency / f0) as usize;
    let signal: Vec<_> = (0..n_step)
        .map(|i| (2f64 * std::f64::consts::PI * i as f64 * f0 / sampling_frequency).sin())
        .collect();

    let tau = sampling_frequency.recip(); // 1/sampling_frequency
                                          // (t,[s])
                                          //let config = complot::Config::new().xaxis(complot::Axis::new().label("Time [s]"));
    let _: complot::Plot = (
        signal
            .iter()
            .enumerate()
            .map(|(i, s)| (i as f64 * tau, vec![*s])),
        complot::complot!("signal.png", xlabel = "Time [s]"),
    )
        .into();

    let welch = Welch::<One>::new(4, 0.5, &signal);
    println!("{welch:}");
    let psd = welch.periogram();

    let sum_sqr = signal.into_iter().map(|x| x * x).sum::<f64>();

    let freq: Vec<_> = (0..psd.len())
        .map(|i| i as f64 * sampling_frequency)
        .collect();
    let n = psd.len();

    let _: complot::LogLog = (
        psd.iter()
            .enumerate()
            .skip(1)
            .map(|(i, s)| (i as f64 * sampling_frequency * 0.5 / n as f64, vec![*s])),
        complot::complot!("psd.png", xlabel = "Frequency [Hz]"),
    )
        .into();
}
