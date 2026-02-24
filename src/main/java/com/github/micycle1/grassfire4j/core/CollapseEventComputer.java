package com.github.micycle1.grassfire4j.core;

import static com.github.micycle1.grassfire4j.core.Core.ccw;
import static com.github.micycle1.grassfire4j.core.Core.cw;
import static com.github.micycle1.grassfire4j.geom.Geom.STOP_EPS;
import static com.github.micycle1.grassfire4j.geom.Geom.dist2;
import static com.github.micycle1.grassfire4j.geom.Geom.getUniqueTimes;
import static com.github.micycle1.grassfire4j.geom.Geom.nearZero;

import java.util.Arrays;

import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.grassfire4j.model.Model.Event;
import com.github.micycle1.grassfire4j.model.Model.Event.EventType;
import com.github.micycle1.grassfire4j.model.Model.KineticTriangle;
import com.github.micycle1.grassfire4j.model.Model.KineticVertex;
import com.github.micycle1.grassfire4j.model.Model.VertexRef;

/**
 * Computes candidate collapse/flip/split events for kinetic triangles.
 */
public class CollapseEventComputer {

	// side indices are only 0, 1, or 2, so can encode them in an int bitmask
	private static final int ALL_SIDES_MASK = 0b111;

	public static double findGt(double[] a, double x) {
		double min = Double.NaN;
		for (int i = 0; i < a.length; i++) {
			double v = a[i];
			if (!Double.isNaN(v) && !nearZero(v - x) && v > x) {
				if (Double.isNaN(min) || v < min) {
					min = v;
				}
			}
		}
		return min;
	}

	public static double findGte(double[] a, double x) {
	    double min = Double.NaN;
	    for (int i = 0; i < a.length; i++) {
	        double v = a[i];
	        if (!Double.isNaN(v) && v >= x) {
	            if (Double.isNaN(min) || v < min) {
	                min = v;
	            }
	        }
	    }
	    return min;
	}

	public static double vertexCrashTime(KineticVertex org, KineticVertex dst, KineticVertex apx, double now) {
		double edgeSpeed = org.wfr.weight;
		Vector2D n = org.ur.w();
		Vector2D s = apx.velocityAt(now);
		double distVE = apx.positionAt(now).subtract(org.positionAt(now)).dot(n);
		double sProj = s.dot(n);
		double denom = edgeSpeed - sProj;
		if (nearZero(denom)) {
			return Double.NaN;
		}
		double tau = distVE / denom;
		if (tau < -STOP_EPS) {
			return Double.NaN;
		}
		return now + Math.max(0.0, tau);
	}

	public static double[] solveQuadratic(double A2, double A1, double A0) {
		if (nearZero(A2) && !nearZero(A1)) {
			return new double[] { -A0 / A1 };
		}
		if (nearZero(A2) && nearZero(A1)) {
			return new double[0];
		}
		double T = -A1 / A2, D = A0 / A2, centre = T * 0.5;
		double under = 0.25 * (T * T) - D;
		if (nearZero(under)) {
			return new double[] { centre };
		}
		if (under < 0) {
			return new double[0];
		}
		double s = Math.sqrt(under);
		double r1 = centre > 0 ? centre + s : centre - s;
		double r2 = r1 != 0 ? D / r1 : 0.0;
		return r1 <= r2 ? new double[] { r1, r2 } : new double[] { r2, r1 };
	}

	public static double[] areaCollapseTimes(VertexRef o, VertexRef d, VertexRef a, double now) {
		Vector2D po = o.positionAt(now), pd = d.positionAt(now), pa = a.positionAt(now);
		Vector2D vo = o.velocityAt(now), vd = d.velocityAt(now), va = a.velocityAt(now);

		double dp10x = pd.getX() - po.getX(), dp10y = pd.getY() - po.getY();
		double dp20x = pa.getX() - po.getX(), dp20y = pa.getY() - po.getY();
		double dv10x = vd.getX() - vo.getX(), dv10y = vd.getY() - vo.getY();
		double dv20x = va.getX() - vo.getX(), dv20y = va.getY() - vo.getY();

		double A2 = dv10x * dv20y - dv10y * dv20x;
		double A1 = (dp10x * dv20y - dp10y * dv20x) + (dv10x * dp20y - dv10y * dp20x);
		double A0 = dp10x * dp20y - dp10y * dp20x;

		double[] out = new double[2];
		int count = 0;
		for (double tau : solveQuadratic(A2, A1, A0)) {
			if (tau >= -STOP_EPS) {
				out[count++] = now + Math.max(0.0, tau);
			}
		}
		if (count == 0) {
			return new double[0];
		}
		if (count == 1) {
			return new double[] { out[0] };
		}
		return out[0] <= out[1] ? out : new double[] { out[1], out[0] };
	}

	public static double collapseTimeEdge(VertexRef v1, VertexRef v2, double now) {
		Vector2D dv = v1.velocityAt(now).subtract(v2.velocityAt(now));
		double denom = dv.dot(dv);
		if (nearZero(denom)) {
			return Double.NaN;
		}
		double tau = dv.dot(v2.positionAt(now).subtract(v1.positionAt(now))) / denom;
		return tau < -STOP_EPS ? Double.NaN : now + Math.max(0.0, tau);
	}

	public static Event compute(KineticTriangle tri, double now, boolean strictGt) {
		if ((tri.stopsAt != null) || !tri.isFinite()) {
			return null; // Logic omitted for infinite to save space. We use bounded CDT.
		}

		int type = tri.getType();
		Event e = null;
		if (type == 0) {
			e = finite0(tri, now, strictGt);
		} else if (type == 1) {
			e = finite1(tri, now, strictGt);
		} else if (type == 2) {
			e = finite2(tri, now, strictGt);
		} else if (type == 3) {
			e = finite3(tri, now, strictGt);
		}

		if (e != null) {
			tri.event = e;
		}
		return e;
	}

	private static void sideD2(KineticVertex o, KineticVertex d, KineticVertex a, double t, double[] out) {
		Vector2D po = o.positionAt(t);
		Vector2D pd = d.positionAt(t);
		Vector2D pa = a.positionAt(t);
		out[0] = dist2(pd, pa);
		out[1] = dist2(pa, po);
		out[2] = dist2(po, pd);
	}

	private static int sideMask(int side) { return 1 << side; }

	private static Event edgeEvt(KineticTriangle t, double time, int sideMask) { return new Event(time, t, sideMask, EventType.EDGE, t.getType()); }
	private static Event flipEvt(KineticTriangle t, double time, int side) { return new Event(time, t, sideMask(side), EventType.FLIP, t.getType()); }
	private static Event splitEvt(KineticTriangle t, double time, int side) { return new Event(time, t, sideMask(side), EventType.SPLIT, t.getType()); }

	private static Event finite0(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		double[] d2 = new double[3];
		double[] tArea = areaCollapseTimes(o, d, a, now);
		for (double t : tArea) {
			if (nearZero(Math.abs(t - now))) {
				sideD2(o, d, a, now, d2);
				int zerosMask = 0;
				for (int i=0; i<3; i++) {
					if (nearZero(Math.sqrt(d2[i]))) {
						zerosMask |= sideMask(i);
					}
				}
				int zerosCount = Integer.bitCount(zerosMask);
				if (zerosCount == 1) {
					return edgeEvt(tri, now, zerosMask);
				}
				if (zerosCount == 3) {
					throw new IllegalStateException("0-triangle collapsing to point at now");
				}
				int mx = d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2);
				return flipEvt(tri, now, mx);
			}
		}

		double[] eTimes = new double[3];
		Arrays.fill(eTimes, Double.NaN);
		for (int i=0; i<3; i++) {
			VertexRef p = i==0?d:(i==1?a:o), q = i==0?a:(i==1?o:d);
			double t = collapseTimeEdge(p, q, now);
			if (!Double.isNaN(t) && nearZero(Math.sqrt(p.distance2At(q, t)))) {
				eTimes[i] = t;
			}
		}

		double te = strictGt ? findGt(eTimes, now) : findGte(eTimes, now);
		double ta = strictGt ? findGt(tArea, now) : findGte(tArea, now);
		if (Double.isNaN(te) && Double.isNaN(ta)) {
			return null;
		}
		if (!Double.isNaN(te) && !Double.isNaN(ta) && nearZero(Math.abs(te - ta))) {
			// Python tie-break: try classify as edge by relative-min; else flip.
			sideD2(o, d, a, te, d2);
			double m = Math.min(d2[0], Math.min(d2[1], d2[2]));
			int relMinMask = 0;
			for (int i=0; i<3; i++) {
				if (nearZero(d2[i] - m)) {
					relMinMask |= sideMask(i);
				}
			}
			int relMinCount = Integer.bitCount(relMinMask);
			if (relMinCount == 3) {
				return edgeEvt(tri, te, ALL_SIDES_MASK);
			}
			if (relMinCount == 1) {
				return edgeEvt(tri, te, relMinMask);
			}
			// otherwise fall back to flip using area time (same)
			sideD2(o, d, a, ta, d2);
			int mx = d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2);
			return flipEvt(tri, ta, mx);
		}

		if (!Double.isNaN(ta) && (Double.isNaN(te) || ta < te)) {
			sideD2(o, d, a, ta, d2);
			return flipEvt(tri, ta, d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2));
		}
		sideD2(o, d, a, te, d2);
		int zerosMask = 0;
		for (int i=0; i<3; i++) {
			if (nearZero(d2[i])) {
				zerosMask |= sideMask(i);
			}
		}
		int zerosCount = Integer.bitCount(zerosMask);
		if (zerosCount == 3) {
			return edgeEvt(tri, te, ALL_SIDES_MASK);
		}
		if (zerosCount == 1) {
			return edgeEvt(tri, te, zerosMask);
		}
		throw new IllegalStateException("0-triangle: edge-collapse time computed but not exactly 1 or 3 sides collapsed");
	}

	private static Event finite1(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		double[] d2 = new double[3];
		int wfSide = tri.indexOfNeighbour(null);
		KineticVertex ow = (KineticVertex)tri.vertices[ccw(wfSide)], dw = (KineticVertex)tri.vertices[cw(wfSide)], aw = (KineticVertex)tri.vertices[wfSide];

		double tc = vertexCrashTime(ow, dw, aw, now);

		// if apex hits wavefront edge "now", classify immediately (edge / split / flip).
		if (!Double.isNaN(tc) && nearZero(Math.abs(tc - now))) {
			sideD2(o, d, a, now, d2);
			int zeroLenMask = 0;
			for (int i = 0; i < 3; i++) {
				if (nearZero(Math.sqrt(d2[i]))) {
					zeroLenMask |= sideMask(i);
				}
			}
			if (Integer.bitCount(zeroLenMask) == 1) {
				return edgeEvt(tri, now, zeroLenMask);
			}
			// else: choose longest side by length; split if it's the wavefront side, else flip.
			double l0 = Math.sqrt(d2[0]), l1 = Math.sqrt(d2[1]), l2 = Math.sqrt(d2[2]);
			int longest = (l0 > l1) ? (l0 > l2 ? 0 : 2) : (l1 > l2 ? 1 : 2);
			return (longest == wfSide) ? splitEvt(tri, now, longest) : flipEvt(tri, now, longest);
		}
		double ta = strictGt ? findGt(areaCollapseTimes(o, d, a, now), now) : findGte(areaCollapseTimes(o, d, a, now), now);
		double te = strictGt ? findGt(new double[] { collapseTimeEdge(ow, dw, now) }, now) : findGte(new double[] { collapseTimeEdge(ow, dw, now) }, now);

		if (Double.isNaN(te) && Double.isNaN(tc)) {
			if (Double.isNaN(ta)) {
				return null;
			}
			if (nearZero(ta - now)) {
				return splitEvt(tri, now, wfSide);
			}
			sideD2(o, d, a, ta, d2);
			for (int i=0; i<3; i++) {
				if (tri.neighbours[i] == null) {
					d2[i] = -1.0;
				}
			}
			return flipEvt(tri, ta, d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2));
		}
		if (Double.isNaN(te) && !Double.isNaN(tc)) {
			double t = (!Double.isNaN(ta) && ta < tc) ? ta : tc;
			sideD2(o, d, a, t, d2);
			double mx = Math.max(d2[0], Math.max(d2[1], d2[2]));
			int longMask = 0;
			for (int i=0; i<3; i++) {
				if (nearZero(Math.sqrt(d2[i]) - Math.sqrt(mx))) {
					longMask |= sideMask(i);
				}
			}
			if (longMask == sideMask(wfSide)) {
				return splitEvt(tri, tc, wfSide);
			}
			int zCt = 0, zIdx = -1;
			for (int i=0; i<3; i++) {
				if (nearZero(Math.sqrt(d2[i]))) { zCt++; zIdx = i; }
			}
			if (zCt == 1) {
				return edgeEvt(tri, t, sideMask(zIdx));
			}
			return flipEvt(tri, t, d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2));
		}
		if (!Double.isNaN(te) && Double.isNaN(tc)) {
			return edgeEvt(tri, te, sideMask(wfSide));
		}
		if (te <= tc) {
			sideD2(o, d, a, te, d2);
			return edgeEvt(tri, te, sideMask(d2[0] < d2[1] ? (d2[0] < d2[2] ? 0 : 2) : (d2[1] < d2[2] ? 1 : 2)));
		}
		sideD2(o, d, a, tc, d2);
		int zCt = 0, zIdx = -1;
		for (int i=0; i<3; i++) {
			if (nearZero(Math.sqrt(d2[i]))) { zCt++; zIdx = i; }
		}
		if (zCt == 1) {
			return edgeEvt(tri, tc, sideMask(zIdx));
		}
		if (zCt == 3) {
			return edgeEvt(tri, tc, ALL_SIDES_MASK);
		}
		int mx = d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2);
		return tri.neighbours[mx] == null ? splitEvt(tri, tc, mx) : flipEvt(tri, tc, mx);
	}

	private static Event finite2(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		double[] l = new double[3];
		double[] tms = new double[] { Double.NaN, Double.NaN, Double.NaN };
		for (int i=0; i<3; i++) {
			if (tri.neighbours[i] == null) {
				tms[i] = collapseTimeEdge(i==0?d:(i==1?a:o), i==0?a:(i==1?o:d), now);
			}
		}
		double t = strictGt ? findGt(getUniqueTimes(tms), now) : findGte(getUniqueTimes(tms), now);
		if (Double.isNaN(t)) {
			t = strictGt ? findGt(areaCollapseTimes(o, d, a, now), now) : findGte(areaCollapseTimes(o, d, a, now), now);
		}
		if (Double.isNaN(t)) {
			return null;
		}
		sideD2(o, d, a, t, l);
		for(int i=0; i<3; i++) {
			l[i] = Math.sqrt(l[i]);
		}
		double mn = Math.min(l[0], Math.min(l[1], l[2]));
		int zCt = 0, zIdx = -1;
		for (int i=0; i<3; i++) {
			if (nearZero(l[i] - mn)) { zCt++; zIdx = i; }
		}
		if (zCt == 3) {
			return edgeEvt(tri, t, ALL_SIDES_MASK);
		}
		if (zCt == 1) {
			return edgeEvt(tri, t, sideMask(zIdx));
		}
		if (zCt == 2) {
			throw new IllegalStateException("2-triangle: impossible configuration (two sides are equal-min)");
		}
		return null; // zCt == 0: no edge event classified
	}

	private static Event finite3(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		// Python: compute edge times per side and which sides are actually ~0 at their closest approach.
		double[] tE = new double[3];
		double[] d2AtTe = new double[3];
		for (int side=0; side<3; side++) {
			VertexRef p = side == 0 ? d : (side == 1 ? a : o);
			VertexRef q = side == 0 ? a : (side == 1 ? o : d);
			tE[side] = collapseTimeEdge(p, q, now);
			d2AtTe[side] = Double.isNaN(tE[side]) ? Double.POSITIVE_INFINITY : p.distance2At(q, tE[side]);
		}

		double tEdge = strictGt ? findGt(tE, now) : findGte(tE, now);
		double tArea = strictGt ? findGt(areaCollapseTimes(o, d, a, now), now) : findGte(areaCollapseTimes(o, d, a, now), now);

		if (!Double.isNaN(tEdge)) {
			int idxMask = 0;
			for (int i=0; i<3; i++) {
				if (nearZero(Math.sqrt(d2AtTe[i]))) {
					idxMask |= sideMask(i);
				}
			}
			int sideMask = idxMask == 0 ? ALL_SIDES_MASK : idxMask;
			int sideCount = Integer.bitCount(sideMask);
			if (sideCount == 2 || sideCount == 0) {
				sideMask = ALL_SIDES_MASK;
			}
			return edgeEvt(tri, tEdge, sideMask);
		}
		if (!Double.isNaN(tArea)) {
			return edgeEvt(tri, tArea, ALL_SIDES_MASK);
		}
		return null;
	}

	public static Event computeNewEdgeCollapse(KineticTriangle tri, double time) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		double[] l = new double[3];
		sideD2(o, d, a, time, l);
		for(int i=0; i<3; i++) {
			l[i] = Math.sqrt(l[i]);
		}
		double mn = Math.min(l[0], Math.min(l[1], l[2]));
		int mask = 0;
		for (int i=0; i<3; i++) {
			if (nearZero(l[i] - mn)) {
				mask |= sideMask(i);
			}
		}
		return edgeEvt(tri, time, mask);
	}
}
