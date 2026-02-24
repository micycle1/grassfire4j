package grassfire2.geom;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class Geom {

	public static record Vec2(double x, double y) {}

	public static final double STOP_EPS = 1e-9;
	public static final double COS_179_999999 = Math.cos(Math.toRadians(179.999999));

	public static double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
	public static Vec2 sub(Vec2 a, Vec2 b) { return new Vec2(a.x - b.x, a.y - b.y); }
	public static Vec2 add(Vec2 a, Vec2 b) { return new Vec2(a.x + b.x, a.y + b.y); }
	public static Vec2 mul(Vec2 a, double s) { return new Vec2(a.x * s, a.y * s); }
	public static double norm2(Vec2 a) { return a.x * a.x + a.y * a.y; }
	public static double norm(Vec2 a) { return Math.hypot(a.x, a.y); }
	public static double dist(Vec2 start, Vec2 end) { return norm(sub(end, start)); }
	public static double dist2(Vec2 start, Vec2 end) { return norm2(sub(end, start)); }

	public static boolean nearZero(double val) {
		return Math.abs(val) <= 1e-10;
	}

	public static record Line2(Vec2 w, double b) {
		public Line2 {
			double n = Math.hypot(w.x(), w.y());
			if (n == 0.0) {
				throw new IllegalArgumentException("Degenerate line with zero normal");
			}
			w = new Vec2(w.x() / n, w.y() / n);
			b = b / n;
		}

		public static Line2 fromPoints(Vec2 start, Vec2 end) {
			double a = start.y() - end.y();
			double bb = end.x() - start.x();
			double c = -start.x() * a - start.y() * bb;
			return new Line2(new Vec2(a, bb), c);
		}

		public Vec2 intersectAtTimeWeighted(Line2 other, double t, double wSelf, double wOther) {
			double a1 = this.w.x(), b1 = this.w.y();
			double a2 = other.w.x(), b2 = other.w.y();
			double c1 = this.b - wSelf * t;
			double c2 = other.b - wOther * t;
			double denom = a1 * b2 - a2 * b1;
			if (nearZero(denom)) {
				return null;
			}
			return new Vec2((b1 * c2 - b2 * c1) / denom, (a2 * c1 - a1 * c2) / denom);
		}
	}

	public static class WaveFront {
		public final Vec2 start, end;
		public final Line2 line;
		public final double weight;
		public final Object data;

		public WaveFront(Vec2 start, Vec2 end, Line2 line, double weight, Object data) {
			this.start = start;
			this.end = end;
			this.line = line != null ? line : Line2.fromPoints(start, end);
			this.weight = weight;
			this.data = data;
		}
	}

	public static Vec2 getBisector(WaveFront left, WaveFront right) {
		Line2 l1 = left.line, l2 = right.line;
		Vec2 p0 = l1.intersectAtTimeWeighted(l2, 0.0, left.weight, right.weight);
		Vec2 p1 = l1.intersectAtTimeWeighted(l2, 1.0, left.weight, right.weight);
		if (p0 == null || p1 == null) {
			double a1 = l1.w.x(), b1 = l1.w.y(), c1 = l1.b;
			double a2 = l2.w.x(), b2 = l2.w.y(), c2 = l2.b;
			double wl = left.weight, wr = right.weight;
			double x1 = a1 * c2 - a2 * c1, x2 = b1 * c2 - b2 * c1;
			if (nearZero(x1) && nearZero(x2)) {
				return new Vec2(0.5 * (wl * a1 + wr * a2), 0.5 * (wl * b1 + wr * b2));
			}
			return new Vec2(wl * a1 + wr * a2, wl * b1 + wr * b2);
		}
		return sub(p1, p0);
	}

	/** Intersection point of the two wavefront support lines at absolute time t (weighted). */
	public static Vec2 intersectionAtTimeWeighted(WaveFront left, WaveFront right, double t) {
		if (left == null || right == null) {
			return null;
		}
		return left.line.intersectAtTimeWeighted(right.line, t, left.weight, right.weight);
	}

	public static List<Double> getUniqueTimes(List<Double> times) {
		if (times == null || times.isEmpty()) {
			return new ArrayList<>();
		}
		List<Double> sorted = new ArrayList<>();
		for (Double d : times) {
			if (d != null) {
				sorted.add(d);
			}
		}
		if (sorted.isEmpty()) {
			return new ArrayList<>();
		}
		Collections.sort(sorted);
		List<Double> out = new ArrayList<>();
		double first = sorted.get(0);
		out.add(first);
		for (double val : sorted) {
			double diff = Math.abs(val - first);
			if (diff <= Math.abs(0.0 * (first + val) * 0.5) || diff <= 1e-7) {
				continue;
			}
			out.add(val);
			first = val;
		}
		return out;
	}
}
