package com.github.micycle1.grassfire4j.core;

import static com.github.micycle1.grassfire4j.core.Core.ccw;
import static com.github.micycle1.grassfire4j.core.Core.cw;
import static com.github.micycle1.grassfire4j.geom.Geom.STOP_EPS;
import static com.github.micycle1.grassfire4j.geom.Geom.dot;
import static com.github.micycle1.grassfire4j.geom.Geom.getUniqueTimes;
import static com.github.micycle1.grassfire4j.geom.Geom.nearZero;
import static com.github.micycle1.grassfire4j.geom.Geom.sub;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.github.micycle1.grassfire4j.geom.Geom.Vec2;
import com.github.micycle1.grassfire4j.model.Model.Event;
import com.github.micycle1.grassfire4j.model.Model.Event.EventType;
import com.github.micycle1.grassfire4j.model.Model.KineticTriangle;
import com.github.micycle1.grassfire4j.model.Model.KineticVertex;
import com.github.micycle1.grassfire4j.model.Model.VertexRef;

/**
 * Computes candidate collapse/flip/split events for kinetic triangles.
 */
public class CollapseEventComputer {

	public static Double findGt(List<Double> a, double x) {
		return a.stream().filter(v -> v != null && !nearZero(v - x) && v > x).min(Double::compare).orElse(null);
	}

	public static Double findGte(List<Double> a, double x) {
		return a.stream().filter(v -> v != null && v >= x).min(Double::compare).orElse(null);
	}

	public static Double vertexCrashTime(KineticVertex org, KineticVertex dst, KineticVertex apx, double now) {
		double edgeSpeed = org.wfr.weight;
		Vec2 n = org.ur.w();
		Vec2 s = apx.velocityAt(now);
		double distVE = dot(sub(apx.positionAt(now), org.positionAt(now)), n);
		double sProj = dot(s, n);
		double denom = edgeSpeed - sProj;
		if (nearZero(denom)) {
			return null;
		}
		double tau = distVE / denom;
		if (tau < -STOP_EPS) {
			return null;
		}
		return now + Math.max(0.0, tau);
	}

	public static List<Double> solveQuadratic(double A2, double A1, double A0) {
		if (nearZero(A2) && !nearZero(A1)) {
			return List.of(-A0 / A1);
		}
		if (nearZero(A2) && nearZero(A1)) {
			return List.of();
		}
		double T = -A1 / A2, D = A0 / A2, centre = T * 0.5;
		double under = 0.25 * (T * T) - D;
		if (nearZero(under)) {
			return List.of(centre);
		}
		if (under < 0) {
			return List.of();
		}
		double s = Math.sqrt(under);
		double r1 = centre > 0 ? centre + s : centre - s;
		double r2 = r1 != 0 ? D / r1 : 0.0;
		List<Double> res = Arrays.asList(r1, r2);
		Collections.sort(res);
		return res;
	}

	public static List<Double> areaCollapseTimes(VertexRef o, VertexRef d, VertexRef a, double now) {
		Vec2 po = o.positionAt(now), pd = d.positionAt(now), pa = a.positionAt(now);
		Vec2 vo = o.velocityAt(now), vd = d.velocityAt(now), va = a.velocityAt(now);

		double dp10x = pd.x() - po.x(), dp10y = pd.y() - po.y();
		double dp20x = pa.x() - po.x(), dp20y = pa.y() - po.y();
		double dv10x = vd.x() - vo.x(), dv10y = vd.y() - vo.y();
		double dv20x = va.x() - vo.x(), dv20y = va.y() - vo.y();

		double A2 = dv10x * dv20y - dv10y * dv20x;
		double A1 = (dp10x * dv20y - dp10y * dv20x) + (dv10x * dp20y - dv10y * dp20x);
		double A0 = dp10x * dp20y - dp10y * dp20x;

		List<Double> out = new ArrayList<>();
		for (double tau : solveQuadratic(A2, A1, A0)) {
			if (tau >= -STOP_EPS) {
				out.add(now + Math.max(0.0, tau));
			}
		}
		Collections.sort(out);
		return out;
	}

	public static Double collapseTimeEdge(VertexRef v1, VertexRef v2, double now) {
		Vec2 dv = sub(v1.velocityAt(now), v2.velocityAt(now));
		double denom = dot(dv, dv);
		if (nearZero(denom)) {
			return null;
		}
		double tau = dot(dv, sub(v2.positionAt(now), v1.positionAt(now))) / denom;
		return tau < -STOP_EPS ? null : now + Math.max(0.0, tau);
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

	private static double[] sideD2(KineticVertex o, KineticVertex d, KineticVertex a, double t) {
		return new double[]{ d.distance2At(a, t), a.distance2At(o, t), o.distance2At(d, t) };
	}

	private static Event edgeEvt(KineticTriangle t, double time, List<Integer> sides) { return new Event(time, t, sides, EventType.EDGE, t.getType()); }
	private static Event flipEvt(KineticTriangle t, double time, int side) { return new Event(time, t, List.of(side), EventType.FLIP, t.getType()); }
	private static Event splitEvt(KineticTriangle t, double time, int side) { return new Event(time, t, List.of(side), EventType.SPLIT, t.getType()); }

	private static Event finite0(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		List<Double> tArea = areaCollapseTimes(o, d, a, now);
		for (double t : tArea) {
			if (nearZero(Math.abs(t - now))) {
				double[] d2 = sideD2(o, d, a, now);
				List<Integer> zeros = new ArrayList<>();
				for (int i=0; i<3; i++) {
					if (nearZero(Math.sqrt(d2[i]))) {
						zeros.add(i);
					}
				}
				if (zeros.size() == 1) {
					return edgeEvt(tri, now, zeros);
				}
				if (zeros.size() == 3) {
					throw new IllegalStateException("0-triangle collapsing to point at now");
				}
				int mx = d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2);
				return flipEvt(tri, now, mx);
			}
		}

		List<Double> eTimes = new ArrayList<>();
		for (int i=0; i<3; i++) {
			VertexRef p = i==0?d:(i==1?a:o), q = i==0?a:(i==1?o:d);
			Double t = collapseTimeEdge(p, q, now);
			if (t != null && nearZero(Math.sqrt(p.distance2At(q, t)))) {
				eTimes.add(t);
			}
		}

		Double te = strictGt ? findGt(eTimes, now) : findGte(eTimes, now);
		Double ta = strictGt ? findGt(tArea, now) : findGte(tArea, now);
		if (te == null && ta == null) {
			return null;
		}
		if (te != null && ta != null && nearZero(Math.abs(te - ta))) {
			// Python tie-break: try classify as edge by relative-min; else flip.
			double[] d2t = sideD2(o, d, a, te);
			double m = Math.min(d2t[0], Math.min(d2t[1], d2t[2]));
			List<Integer> relMin = new ArrayList<>();
			for (int i=0; i<3; i++) {
				if (nearZero(d2t[i] - m)) {
					relMin.add(i);
				}
			}
			if (relMin.size() == 3) {
				return edgeEvt(tri, te, List.of(0,1,2));
			}
			if (relMin.size() == 1) {
				return edgeEvt(tri, te, List.of(relMin.get(0)));
			}
			// otherwise fall back to flip using area time (same)
			double[] d2a = sideD2(o, d, a, ta);
			int mx = d2a[0] > d2a[1] ? (d2a[0] > d2a[2] ? 0 : 2) : (d2a[1] > d2a[2] ? 1 : 2);
			return flipEvt(tri, ta, mx);
		}

		if (ta != null && (te == null || ta < te)) {
			double[] d2 = sideD2(o, d, a, ta);
			return flipEvt(tri, ta, d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2));
		}
		double[] d2 = sideD2(o, d, a, te);
		List<Integer> zeros = new ArrayList<>();
		for (int i=0; i<3; i++) {
			if (nearZero(d2[i])) {
				zeros.add(i);
			}
		}
		if (zeros.size() == 3) {
			return edgeEvt(tri, te, List.of(0,1,2));
		}
		if (zeros.size() == 1) {
			return edgeEvt(tri, te, List.of(zeros.get(0)));
		}
		throw new IllegalStateException("0-triangle: edge-collapse time computed but not exactly 1 or 3 sides collapsed");
	}

	private static Event finite1(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		int wfSide = tri.indexOfNeighbour(null);
		KineticVertex ow = (KineticVertex)tri.vertices[ccw(wfSide)], dw = (KineticVertex)tri.vertices[cw(wfSide)], aw = (KineticVertex)tri.vertices[wfSide];

		Double tc = vertexCrashTime(ow, dw, aw, now);

		// if apex hits wavefront edge "now", classify immediately (edge / split / flip).
		if (tc != null && nearZero(Math.abs(tc - now))) {
			double[] d2Now = sideD2(o, d, a, now);
			List<Integer> zeroLenSides = new ArrayList<>();
			for (int i = 0; i < 3; i++) {
				if (nearZero(Math.sqrt(d2Now[i]))) {
					zeroLenSides.add(i);
				}
			}
			if (zeroLenSides.size() == 1) {
				return edgeEvt(tri, now, List.of(zeroLenSides.get(0)));
			}
			// else: choose longest side by length; split if it's the wavefront side, else flip.
			double l0 = Math.sqrt(d2Now[0]), l1 = Math.sqrt(d2Now[1]), l2 = Math.sqrt(d2Now[2]);
			int longest = (l0 > l1) ? (l0 > l2 ? 0 : 2) : (l1 > l2 ? 1 : 2);
			return (longest == wfSide) ? splitEvt(tri, now, longest) : flipEvt(tri, now, longest);
		}
		Double ta = strictGt ? findGt(areaCollapseTimes(o, d, a, now), now) : findGte(areaCollapseTimes(o, d, a, now), now);
		Double te = strictGt ? findGt(Arrays.asList(collapseTimeEdge(ow, dw, now)), now) : findGte(Arrays.asList(collapseTimeEdge(ow, dw, now)), now);

		if (te == null && tc == null) {
			if (ta == null) {
				return null;
			}
			if (nearZero(ta - now)) {
				return splitEvt(tri, now, wfSide);
			}
			double[] d2 = sideD2(o, d, a, ta);
			for (int i=0; i<3; i++) {
				if (tri.neighbours[i] == null) {
					d2[i] = -1.0;
				}
			}
			return flipEvt(tri, ta, d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2));
		}
		if (te == null && tc != null) {
			double t = (ta != null && ta < tc) ? ta : tc;
			double[] d2 = sideD2(o, d, a, t);
			double mx = Math.max(d2[0], Math.max(d2[1], d2[2]));
			List<Integer> longs = new ArrayList<>();
			for (int i=0; i<3; i++) {
				if (nearZero(Math.sqrt(d2[i]) - Math.sqrt(mx))) {
					longs.add(i);
				}
			}
			if (longs.size() == 1 && longs.get(0) == wfSide) {
				return splitEvt(tri, tc, wfSide);
			}
			int zCt = 0, zIdx = -1;
			for (int i=0; i<3; i++) {
				if (nearZero(Math.sqrt(d2[i]))) { zCt++; zIdx = i; }
			}
			if (zCt == 1) {
				return edgeEvt(tri, t, List.of(zIdx));
			}
			return flipEvt(tri, t, d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2));
		}
		if (te != null && tc == null) {
			return edgeEvt(tri, te, List.of(wfSide));
		}
		if (te <= tc) {
			double[] d2 = sideD2(o, d, a, te);
			return edgeEvt(tri, te, List.of(d2[0] < d2[1] ? (d2[0] < d2[2] ? 0 : 2) : (d2[1] < d2[2] ? 1 : 2)));
		}
		double[] d2 = sideD2(o, d, a, tc);
		int zCt = 0, zIdx = -1;
		for (int i=0; i<3; i++) {
			if (nearZero(Math.sqrt(d2[i]))) { zCt++; zIdx = i; }
		}
		if (zCt == 1) {
			return edgeEvt(tri, tc, List.of(zIdx));
		}
		if (zCt == 3) {
			return edgeEvt(tri, tc, List.of(0,1,2));
		}
		int mx = d2[0] > d2[1] ? (d2[0] > d2[2] ? 0 : 2) : (d2[1] > d2[2] ? 1 : 2);
		return tri.neighbours[mx] == null ? splitEvt(tri, tc, mx) : flipEvt(tri, tc, mx);
	}

	private static Event finite2(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		List<Double> tms = new ArrayList<>();
		for (int i=0; i<3; i++) {
			if (tri.neighbours[i] == null) {
				tms.add(collapseTimeEdge(i==0?d:(i==1?a:o), i==0?a:(i==1?o:d), now));
			}
		}
		Double t = strictGt ? findGt(getUniqueTimes(tms), now) : findGte(getUniqueTimes(tms), now);
		if (t == null) {
			t = strictGt ? findGt(areaCollapseTimes(o, d, a, now), now) : findGte(areaCollapseTimes(o, d, a, now), now);
		}
		if (t == null) {
			return null;
		}
		double[] l = sideD2(o, d, a, t);
		for(int i=0; i<3; i++) {
			l[i] = Math.sqrt(l[i]);
		}
		double mn = Math.min(l[0], Math.min(l[1], l[2]));
		int zCt = 0, zIdx = -1;
		for (int i=0; i<3; i++) {
			if (nearZero(l[i] - mn)) { zCt++; zIdx = i; }
		}
		if (zCt == 3) {
			return edgeEvt(tri, t, List.of(0,1,2));
		}
		if (zCt == 1) {
			return edgeEvt(tri, t, List.of(zIdx));
		}
		if (zCt == 2) {
			throw new IllegalStateException("2-triangle: impossible configuration (two sides are equal-min)");
		}
		return null; // zCt == 0: no edge event classified
	}

	private static Event finite3(KineticTriangle tri, double now, boolean strictGt) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		// Python: compute edge times per side and which sides are actually ~0 at their closest approach.
		Double[] tE = new Double[3];
		double[] d2AtTe = new double[3];
		for (int side=0; side<3; side++) {
			VertexRef p = side==0?d:(side==1?a:o);
			VertexRef q = side==0?a:(side==1?o:d);
			tE[side] = collapseTimeEdge(p, q, now);
			d2AtTe[side] = (tE[side] == null) ? Double.POSITIVE_INFINITY : p.distance2At(q, tE[side]);
		}
		List<Double> edgeTimes = new ArrayList<>();
		for (Double t : tE) {
			if (t != null) {
				edgeTimes.add(t);
			}
		}

		Double tEdge = strictGt ? findGt(edgeTimes, now) : findGte(edgeTimes, now);
		Double tArea = strictGt ? findGt(areaCollapseTimes(o, d, a, now), now) : findGte(areaCollapseTimes(o, d, a, now), now);

		if (tEdge != null) {
			List<Integer> idx = new ArrayList<>();
			for (int i=0; i<3; i++) {
				if (nearZero(Math.sqrt(d2AtTe[i]))) {
					idx.add(i);
				}
			}
			List<Integer> sides = idx.isEmpty() ? List.of(0,1,2) : idx;
			if (sides.size() == 2 || sides.isEmpty()) {
				sides = List.of(0,1,2);
			}
			return edgeEvt(tri, tEdge, sides);
		}
		if (tArea != null) {
			return edgeEvt(tri, tArea, List.of(0,1,2));
		}
		return null;
	}

	public static Event computeNewEdgeCollapse(KineticTriangle tri, double time) {
		KineticVertex o = (KineticVertex)tri.vertices[0], d = (KineticVertex)tri.vertices[1], a = (KineticVertex)tri.vertices[2];
		double[] l = sideD2(o, d, a, time);
		for(int i=0; i<3; i++) {
			l[i] = Math.sqrt(l[i]);
		}
		double mn = Math.min(l[0], Math.min(l[1], l[2]));
		List<Integer> sides = new ArrayList<>();
		for (int i=0; i<3; i++) {
			if (nearZero(l[i] - mn)) {
				sides.add(i);
			}
		}
		return edgeEvt(tri, time, sides);
	}
}
