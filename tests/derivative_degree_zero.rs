use rusty_fitpack::{
    splder, splder_uniform, splder_uniform_ext, splev, splev_uniform, Extrapolation,
};

fn monomial(degree: usize) -> (Vec<f64>, Vec<f64>) {
    let mut knots = vec![0.0; degree + 1];
    knots.extend(vec![1.0; degree + 1]);
    let mut coefficients = vec![0.0; degree + 1];
    coefficients[degree] = 1.0;
    (knots, coefficients)
}

#[test]
fn linear_single_interval_has_unit_slope() {
    let (knots, coefficients) = monomial(1);
    for x in [0.0, 0.25, 0.75, 1.0] {
        assert_eq!(splder_uniform(&knots, &coefficients, 1, x, 1), 1.0);
    }
}

#[test]
fn piecewise_linear_uses_current_interval_and_clamps() {
    let knots = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];
    let coefficients = vec![0.0, 1.0, 3.0, 6.0];
    let points = vec![-0.25, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.25];
    let expected = vec![1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0];
    assert_eq!(splder(&knots, &coefficients, 1, &points, 1), expected);
    for (&x, &slope) in points.iter().zip(&expected) {
        assert_eq!(splder_uniform(&knots, &coefficients, 1, x, 1), slope);
    }
}

#[test]
fn linear_derivative_matches_value_only_finite_difference() {
    let (knots, coefficients) = monomial(1);
    let h = 1.0e-5;
    for x in [0.125, 0.5, 0.875] {
        let numerical = (splev_uniform(&knots, &coefficients, 1, x + h)
            - splev_uniform(&knots, &coefficients, 1, x - h))
            / (2.0 * h);
        let analytical = splder_uniform(&knots, &coefficients, 1, x, 1);
        assert!((analytical - numerical).abs() < 1.0e-8);
    }
}

#[test]
fn piecewise_linear_matches_general_value_finite_difference() {
    let knots = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];
    let coefficients = vec![0.0, 1.0, 3.0, 6.0];
    let h = 1.0e-5;
    for x in [0.25, 0.75, 1.25, 1.75, 2.25, 2.75] {
        let values = splev(
            knots.clone(),
            coefficients.clone(),
            1,
            vec![x - h, x + h],
            0,
        );
        let numerical = (values[1] - values[0]) / (2.0 * h);
        let analytical = splder_uniform(&knots, &coefficients, 1, x, 1);
        assert!((analytical - numerical).abs() < 1.0e-8);
    }
}

#[test]
fn highest_derivatives_of_monomials_are_factorials() {
    let mut factorial = 1.0;
    for degree in 1..=5 {
        factorial *= degree as f64;
        let (knots, coefficients) = monomial(degree);
        for x in [-0.25, 0.0, 0.125, 0.75, 1.0, 1.25] {
            assert_eq!(
                splder_uniform(&knots, &coefficients, degree, x, degree),
                factorial
            );
        }
    }
}

#[test]
fn cubic_third_derivative_uses_current_interval() {
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0];
    let coefficients = vec![0.0, 1.0, 4.0, 10.0, 20.0, 35.0];
    let points = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
    let expected = vec![-1.5, -1.5, 7.5, 7.5, 51.0, 51.0, 51.0];
    assert_eq!(splder(&knots, &coefficients, 3, &points, 3), expected);
    for (&x, &value) in points.iter().zip(&expected) {
        assert_eq!(splder_uniform(&knots, &coefficients, 3, x, 3), value);
    }
}

#[test]
fn ordinary_cubic_first_and_second_derivatives_remain_exact() {
    let (knots, coefficients) = monomial(3);
    for x in [0.0, 0.125, 0.5, 0.75, 1.0] {
        assert_eq!(splder_uniform(&knots, &coefficients, 3, x, 1), 3.0 * x * x);
        assert_eq!(splder_uniform(&knots, &coefficients, 3, x, 2), 6.0 * x);
    }
    assert_eq!(splder_uniform(&knots, &coefficients, 3, -0.5, 1), 0.0);
    assert_eq!(splder_uniform(&knots, &coefficients, 3, 1.5, 1), 3.0);
}

#[test]
fn explicit_extrapolation_policies_keep_their_contract() {
    let (knots, coefficients) = monomial(1);
    for mode in [
        Extrapolation::Clamp,
        Extrapolation::Zero,
        Extrapolation::SmoothDecay(1.0),
    ] {
        assert_eq!(
            splder_uniform_ext(&knots, &coefficients, 1, 0.75, 1, mode),
            1.0
        );
    }
    for mode in [Extrapolation::Clamp, Extrapolation::Zero] {
        assert_eq!(
            splder_uniform_ext(&knots, &coefficients, 1, -0.5, 1, mode),
            0.0
        );
        assert_eq!(
            splder_uniform_ext(&knots, &coefficients, 1, 1.5, 1, mode),
            0.0
        );
    }
}
