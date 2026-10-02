use rusty_fitpack::{
    splder, splder_uniform, splder_uniform_ext, splev, splev_uniform, splev_uniform_ext, splrep,
    Extrapolation,
};

fn fitted(degree: usize, origin: f64, step: f64) -> (Vec<f64>, Vec<f64>, usize) {
    let x: Vec<f64> = (0..12).map(|i| origin + step * i as f64).collect();
    let y: Vec<f64> = (0..12).map(|i| (-(i as f64).powi(2) / 8.0).exp()).collect();
    splrep(
        x,
        y,
        None,
        None,
        None,
        Some(degree),
        None,
        None,
        None,
        None,
        None,
        None,
    )
}

fn clamped_cubic() -> (Vec<f64>, Vec<f64>) {
    (
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0],
        vec![0.0, 1.0, 4.0, 10.0, 20.0, 35.0],
    )
}

#[test]
fn fitted_linear_values_use_current_interval() {
    let (t, c, k) = splrep(
        vec![0.0, 1.0, 2.0, 3.0],
        vec![0.0, 1.0, 3.0, 6.0],
        None,
        None,
        None,
        Some(1),
        None,
        None,
        None,
        None,
        None,
        None,
    );
    for (x, expected) in [(0.5, 0.5), (1.5, 2.0), (2.5, 4.5)] {
        assert_eq!(splev_uniform(&t, &c, k, x), expected);
        assert_eq!(splder_uniform(&t, &c, k, x, 0), expected);
    }
}

#[test]
fn piecewise_linear_uniform_value_matches_derivative() {
    let t = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];
    let c = vec![0.0, 1.0, 3.0, 6.0];
    let h = 1.0e-5;
    for (x, expected) in [(0.5, 1.0), (1.5, 2.0), (2.5, 3.0)] {
        let numerical =
            (splev_uniform(&t, &c, 1, x + h) - splev_uniform(&t, &c, 1, x - h)) / (2.0 * h);
        assert_eq!(splder_uniform(&t, &c, 1, x, 1), expected);
        assert!((numerical - expected).abs() < 1.0e-8);
    }
}

#[test]
fn clamped_cubic_matches_independent_rational_values() {
    let (t, c) = clamped_cubic();
    // Exact Cox-de Boor results, independently derived with rational arithmetic.
    for (x, value, slope) in [
        (0.5, 59.0 / 32.0, 69.0 / 16.0),
        (1.5, 231.0 / 32.0, 111.0 / 16.0),
        (2.5, 303.0 / 16.0, 171.0 / 8.0),
    ] {
        assert!((splev_uniform(&t, &c, 3, x) - value).abs() < 1.0e-12);
        assert_eq!(splder_uniform(&t, &c, 3, x, 1), slope);
    }
}

#[test]
fn fitted_cubic_first_interval_derivatives_match_general_evaluator() {
    let (t, c, k) = fitted(3, 0.0, 1.0);
    let points = vec![0.0, 0.125, 0.5, 0.875, 1.0, 1.5];
    for order in 1..=2 {
        let expected = splder(&t, &c, k, &points, order);
        for (&x, &value) in points.iter().zip(&expected) {
            assert_eq!(splder_uniform(&t, &c, k, x, order), value);
        }
    }
}

#[test]
fn fitted_cubic_first_interval_matches_value_only_finite_difference() {
    let (t, c, k) = fitted(3, 0.0, 1.0);
    let h = 1.0e-5;
    for x in [0.125, 0.5, 0.875] {
        let numerical =
            (splev_uniform(&t, &c, k, x + h) - splev_uniform(&t, &c, k, x - h)) / (2.0 * h);
        assert!((splder_uniform(&t, &c, k, x, 1) - numerical).abs() < 1.0e-8);
    }
}

#[test]
fn degrees_one_through_five_match_general_evaluators_on_shifted_grids() {
    for degree in 1..=5 {
        for (origin, step) in [(0.0, 1.0), (0.37, 0.02), (-3.0, 0.25)] {
            let (t, c, k) = fitted(degree, origin, step);
            let points: Vec<f64> = (0..=88).map(|i| origin + step * i as f64 / 8.0).collect();
            let expected = splev(t.clone(), c.clone(), k, points.clone(), 3);
            for (&x, &value) in points.iter().zip(&expected) {
                assert_eq!(splev_uniform(&t, &c, k, x), value);
                assert_eq!(splder_uniform(&t, &c, k, x, 0), value);
            }
            for order in 1..=k {
                let expected = splder(&t, &c, k, &points, order);
                for (&x, &value) in points.iter().zip(&expected) {
                    assert_eq!(splder_uniform(&t, &c, k, x, order), value);
                }
            }
        }
    }
}

#[test]
fn values_and_derivatives_match_general_evaluators_around_knots() {
    let (t, c) = clamped_cubic();
    let points = vec![
        0.0,
        1.0 - 1.0e-12,
        1.0,
        1.0 + 1.0e-12,
        2.0 - 1.0e-12,
        2.0,
        2.0 + 1.0e-12,
        3.0,
    ];
    let expected = splev(t.clone(), c.clone(), 3, points.clone(), 3);
    for (&x, &value) in points.iter().zip(&expected) {
        assert_eq!(splev_uniform(&t, &c, 3, x), value);
    }
    for order in 0..=3 {
        let expected = splder(&t, &c, 3, &points, order);
        for (&x, &value) in points.iter().zip(&expected) {
            assert_eq!(splder_uniform(&t, &c, 3, x, order), value);
        }
    }
}

#[test]
fn corrected_intervals_preserve_extrapolation_contracts() {
    let (t, c) = clamped_cubic();
    for mode in [
        Extrapolation::Clamp,
        Extrapolation::Zero,
        Extrapolation::SmoothDecay(1.0),
    ] {
        for x in [0.125, 1.5, 2.5, 3.0] {
            assert_eq!(
                splev_uniform_ext(&t, &c, 3, x, mode),
                splev_uniform(&t, &c, 3, x)
            );
            for order in 0..=3 {
                assert_eq!(
                    splder_uniform_ext(&t, &c, 3, x, order, mode),
                    splder_uniform(&t, &c, 3, x, order)
                );
            }
        }
        assert_eq!(splder_uniform_ext(&t, &c, 3, -0.5, 1, mode), 0.0);
    }
    assert_eq!(
        splev_uniform(&t, &c, 3, -0.5),
        splev_uniform(&t, &c, 3, 0.0)
    );
    assert_eq!(splev_uniform(&t, &c, 3, 3.5), splev_uniform(&t, &c, 3, 3.0));
    for order in 1..=3 {
        for (outside, endpoint) in [(-0.5, 0.0), (3.5, 3.0)] {
            assert_eq!(
                splder_uniform(&t, &c, 3, outside, order),
                splder_uniform(&t, &c, 3, endpoint, order)
            );
        }
        for mode in [Extrapolation::Clamp, Extrapolation::Zero] {
            assert_eq!(splder_uniform_ext(&t, &c, 3, 3.5, order, mode), 0.0);
        }
        assert_eq!(
            splder_uniform_ext(&t, &c, 3, 4.0, order, Extrapolation::SmoothDecay(1.0)),
            0.0
        );
    }
    assert_eq!(splev_uniform_ext(&t, &c, 3, 3.5, Extrapolation::Zero), 0.0);
    assert_eq!(
        splev_uniform_ext(&t, &c, 3, 4.0, Extrapolation::SmoothDecay(1.0)),
        0.0
    );
}
