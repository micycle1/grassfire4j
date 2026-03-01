package com.github.micycle1.grassfire4j.geom;

import java.util.Arrays;

import org.locationtech.jts.math.Vector2D;

/**
 * Geometric primitives and utility operations used by the kinetic solver.
 * <p>
 * Shared maths lives here so adapter code and solver code can rely on a common,
 * minimal geometry layer.
 */
public final class Geom {

	public static final double STOP_EPS = 1e-9;
	private static final double TIME_UNIQUE_TOL = 1e-9;
	private static final double NEAR_ZERO_EPS = 1e-8;
	private static final double NEAR_ZERO_EPS_SQ = NEAR_ZERO_EPS * NEAR_ZERO_EPS;

	private Geom() {
	}

	public static double dist2(Vector2D start, Vector2D end) {
		return end.subtract(start).lengthSquared();
	}

	public static boolean nearZero(double val) {
		return Math.abs(val) <= NEAR_ZERO_EPS;
	}

	public static boolean nearZeroTime(double val) {
		return Math.abs(val) <= TIME_UNIQUE_TOL;
	}

	public static boolean nearZeroSq(double valSq) {
		return Math.abs(valSq) <= NEAR_ZERO_EPS_SQ;
	}

	public static record Line2(Vector2D w, double b) {
		public Line2 {
			double n = Math.hypot(w.getX(), w.getY());
			if (n == 0.0) {
				throw new IllegalArgumentException("Degenerate line with zero normal");
			}
			w = new Vector2D(w.getX() / n, w.getY() / n);
			b = b / n;
		}

		public static Line2 fromPoints(Vector2D start, Vector2D end) {
			double a = start.getY() - end.getY();
			double bb = end.getX() - start.getX();
			double c = -start.getX() * a - start.getY() * bb;
			return new Line2(new Vector2D(a, bb), c);
		}

		public Vector2D intersectAtTimeWeighted(Line2 other, double t, double wSelf, double wOther) {
			double a1 = this.w.getX(), b1 = this.w.getY();
			double a2 = other.w.getX(), b2 = other.w.getY();
			double c1 = this.b - wSelf * t;
			double c2 = other.b - wOther * t;
			double denom = a1 * b2 - a2 * b1;
			if (nearZero(denom)) {
				return null;
			}
			return new Vector2D((b1 * c2 - b2 * c1) / denom, (a2 * c1 - a1 * c2) / denom);
		}
	}

	public static class WaveFront {
		public final Vector2D start, end;
		public final Line2 line;
		public final double weight;
		public final Object data;

		public WaveFront(Vector2D start, Vector2D end, Line2 line, double weight, Object data) {
			this.start = start;
			this.end = end;
			this.line = line != null ? line : Line2.fromPoints(start, end);
			this.weight = weight;
			this.data = data;
		}
	}

	public static Vector2D getBisector(WaveFront left, WaveFront right) {
		Line2 l1 = left.line, l2 = right.line;
		Vector2D p0 = l1.intersectAtTimeWeighted(l2, 0.0, left.weight, right.weight);
		Vector2D p1 = l1.intersectAtTimeWeighted(l2, 1.0, left.weight, right.weight);
		if (p0 == null || p1 == null) {
			double a1 = l1.w.getX(), b1 = l1.w.getY(), c1 = l1.b;
			double a2 = l2.w.getX(), b2 = l2.w.getY(), c2 = l2.b;
			double wl = left.weight, wr = right.weight;
			double x1 = a1 * c2 - a2 * c1, x2 = b1 * c2 - b2 * c1;
			if (nearZero(x1) && nearZero(x2)) {
				return new Vector2D(0.5 * (wl * a1 + wr * a2), 0.5 * (wl * b1 + wr * b2));
			}
			return new Vector2D(wl * a1 + wr * a2, wl * b1 + wr * b2);
		}
		return p1.subtract(p0);
	}

	/**
	 * Intersection point of the two wavefront support lines at absolute time t
	 * (weighted).
	 */
	public static Vector2D intersectionAtTimeWeighted(WaveFront left, WaveFront right, double t) {
		if (left == null || right == null) {
			return null;
		}
		return left.line.intersectAtTimeWeighted(right.line, t, left.weight, right.weight);
	}

	public static double[] getUniqueTimes(double[] times) {
		if (times == null || times.length == 0) {
			return new double[0];
		}
		double[] sorted = new double[times.length];
		int n = 0;
		for (double d : times) {
			if (!Double.isNaN(d)) {
				sorted[n++] = d;
			}
		}
		if (n == 0) {
			return new double[0];
		}
		Arrays.sort(sorted, 0, n);
		double[] out = new double[n];
		int m = 0;
		double first = sorted[0];
		out[m++] = first;
		for (int i = 1; i < n; i++) {
			double val = sorted[i];
			double diff = Math.abs(val - first);
			// use a combined relative/absolute tolerance
			double relTol = Math.abs((first + val) * 0.5) * TIME_UNIQUE_TOL;
			double absTol = TIME_UNIQUE_TOL;
			if (diff <= Math.max(relTol, absTol)) {
				continue;
			}
			out[m++] = val;
			first = val;
		}
		return m == n ? out : Arrays.copyOf(out, m);
	}
}
