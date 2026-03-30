use std::f64::consts::PI;

use plotters::prelude::*;

fn f(x: f64) -> f64 {
    1.0 / (1.0 + 25.0 * x * x)
}

fn xi(n: u32) -> impl Iterator<Item = f64> {
    (0..=n).map(move |i| 2.0 * i as f64 / n as f64 - 1.0)
}

fn chebyshev_points(n: u32, a: f64, b: f64) -> impl Iterator<Item = f64> {
    (0..=n).rev().map(move |k| {
        0.5 * (b + a) + 0.5 * (b - a) * ((2.0 * k as f64 + 1.0) * PI / (2.0 * (n + 1) as f64)).cos()
    })
}

fn divided_difference(points: &[(f64, f64)]) -> f64 {
    if points.len() == 0 {
        panic!("no points");
    }

    if points.len() == 1 {
        points[0].1
    } else {
        let (pi, pjm) = points.split_first().unwrap();
        let (pm, pik) = points.split_last().unwrap();
        (divided_difference(pjm) - divided_difference(pik)) / (pm.0 - pi.0)
    }
}

fn newtons_polynomial(points: &[(f64, f64)]) -> impl Fn(f64) -> f64 {
    let diffs: Vec<_> = (1..points.len())
        .map(|i| divided_difference(&points[0..=i]))
        .collect();

    move |x| {
        let mut res = points[0].1;
        let mut coeff = 1.0;
        for i in 0..points.len() - 1 {
            coeff *= x - points[i].0;
            res += coeff * diffs[i];
        }
        res
    }
}

const HALF_RESOLUTION: u32 = 1000;

struct Function<'a> {
    pub name: &'a str,
    pub f: Box<dyn Fn(f64) -> f64 + 'a>,
}

fn draw<'a, P, S>(path: &P, plot_name: &S, fns: &[Function<'a>]) -> anyhow::Result<()>
where
    P: AsRef<std::path::Path> + ?Sized,
    S: AsRef<str> + ?Sized,
{
    let root = BitMapBackend::new(&path, (1600, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let points_to_draw = fns.iter().map(|f| {
        (-(HALF_RESOLUTION as i64)..=(HALF_RESOLUTION as i64))
            .map(|x| x as f64 / HALF_RESOLUTION as f64)
            .map(|x| (x, f.f.as_ref()(x)))
    });

    let max_error = points_to_draw
        .clone()
        .flat_map(|it| it.map(|t| t.1))
        .max_by(|x, y| x.partial_cmp(y).expect("a valid float"))
        .unwrap();

    let colors = [RED, BLUE, GREEN].iter().cycle();

    let series = points_to_draw
        .clone()
        .zip(colors.clone())
        .map(|(p, c)| LineSeries::new(p, c));

    let mut chart = ChartBuilder::on(&root)
        .margin(16)
        .caption(plot_name, ("serif", 20).into_font())
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-1.0..1.0, 0.0..max_error.ceil())?;

    chart.configure_mesh().x_desc("x").y_desc("error").draw()?;

    for ((series, name), c) in series.zip(fns.iter().map(|f| f.name)).zip(colors) {
        chart
            .draw_series(series)?
            .label(name)
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], *c));
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

fn main() -> anyhow::Result<()> {
    for n in 3..=10 {
        let xi: Vec<_> = xi(n).map(|x| (x, f(x))).collect();
        let cheb: Vec<_> = chebyshev_points(n, -1.0, 1.0).map(|x| (x, f(x))).collect();

        let p1 = newtons_polynomial(&xi);
        let p2 = newtons_polynomial(&cheb);

        let path = format!("plots/{n:02}.png");

        draw(
            &path,
            &format!("n = {n}"),
            &[
                Function {
                    name: "x_n",
                    f: Box::new(|x| (f(x) - p1(x)).abs()),
                },
                Function {
                    name: "chebyshev",
                    f: Box::new(|x| (f(x) - p2(x)).abs()),
                },
            ],
        )?;
    }

    Ok(())
}
