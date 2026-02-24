package com.github.micycle1.grassfire4j.events;

import static com.github.micycle1.grassfire4j.core.Core.apex;
import static com.github.micycle1.grassfire4j.core.Core.ccw;
import static com.github.micycle1.grassfire4j.core.Core.cw;
import static com.github.micycle1.grassfire4j.core.Core.dest;
import static com.github.micycle1.grassfire4j.core.Core.orig;
import static com.github.micycle1.grassfire4j.geom.Geom.COS_179_999999;
import static com.github.micycle1.grassfire4j.geom.Geom.nearZero;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.List;
import java.util.PriorityQueue;
import java.util.function.IntUnaryOperator;

import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.grassfire4j.core.CollapseEventComputer;
import com.github.micycle1.grassfire4j.core.Core;
import com.github.micycle1.grassfire4j.geom.Geom;
import com.github.micycle1.grassfire4j.geom.Geom.Line2;
import com.github.micycle1.grassfire4j.geom.Geom.WaveFront;
import com.github.micycle1.grassfire4j.model.Model.Event;
import com.github.micycle1.grassfire4j.model.Model.Event.EventType;
import com.github.micycle1.grassfire4j.model.Model.KineticTriangle;
import com.github.micycle1.grassfire4j.model.Model.KineticVertex;
import com.github.micycle1.grassfire4j.model.Model.Skeleton;
import com.github.micycle1.grassfire4j.model.Model.SkeletonNode;
import com.github.micycle1.grassfire4j.model.Model.VertexRef;

/**
 * Event queue and handlers for advancing the kinetic straight-skeleton
 * simulation.
 */
public class Events {

	private static final int LOOP_MAX = 50_000;

	public static class EventQueue {
		private final PriorityQueue<Event> heap = new PriorityQueue<>();
		private long counter = 0;

		public void add(Event item) {
			item.valid = true;
			item.counter = counter++;
			heap.add(item);
		}

		public void discard(Event item) { item.valid = false; }

		public Event pop() {
			while (!heap.isEmpty()) {
				Event item = heap.poll();
				if (item.valid) {
					return item;
				}
			}
			return null;
		}

		public boolean isEmpty() {
			while (!heap.isEmpty() && !heap.peek().valid) {
				heap.poll();
			}
			return heap.isEmpty();
		}
	}

	public static void stopKVertices(List<KineticVertex> V, int step, double now, Skeleton skel, Vector2D pos) {
		SkeletonNode skNode = null;
		boolean isNew = false;
		for (var v : V) {
			if (v.stopsAt != null) {
				skNode = v.stopNode;
			} else if (nearZero(v.startsAt - now)) {
				skNode = v.startNode;
			} else {
				v.stopsAt = now;
			}
		}
		if (skNode == null) {
			if (pos == null) {
				double sx=0, sy=0;
				for (var v : V) { Vector2D p = v.positionAt(now); sx += p.getX(); sy += p.getY(); }
				pos = new Vector2D(sx / V.size(), sy / V.size());
			}
			skNode = new SkeletonNode(pos, step, null);
			isNew = true;
		}
		for (var v : V) { v.stopNode = skNode; v.stopsAt = now; }
		if (isNew) {
			skel.skNodes.add(skNode);
		}
	}

	public static KineticVertex computeNewKVertex(WaveFront wfl, WaveFront wfr, double now, SkeletonNode skNode, int info, boolean internal) {
		KineticVertex kv = new KineticVertex();
		kv.info = info; kv.startsAt = now; kv.startNode = skNode; kv.internal = internal;
		Line2 ul = wfl.line, ur = wfr.line;
		double u1x = ul.w().getX(), u1y = ul.w().getY(), u2x = ur.w().getX(), u2y = ur.w().getY();
		double d = Math.max(-1.0, Math.min(1.0, u1x * u2x + u1y * u2y));
		double dirx = wfl.weight * u1x + wfr.weight * u2x, diry = wfl.weight * u1y + wfr.weight * u2y;

		Vector2D bi, pos0;
		if ((nearZero(dirx) && nearZero(diry)) || nearZero(d + 1.0) || d < COS_179_999999) {
			bi = new Vector2D(0, 0); pos0 = skNode.pos;
		} else {
			double denom = u1x * u2y - u2x * u1y;
			if (nearZero(denom)) {
				double x1 = u1x * ur.b() - u2x * ul.b(), x2 = u1y * ur.b() - u2y * ul.b();
				if (nearZero(x1) && nearZero(x2)) {
					bi = new Vector2D(0.5 * dirx, 0.5 * diry);
					pos0 = new Vector2D(skNode.pos.getX() - bi.getX() * now, skNode.pos.getY() - bi.getY() * now);
				} else { bi = new Vector2D(0, 0); pos0 = skNode.pos; }
			} else {
				Vector2D p0 = ul.intersectAtTimeWeighted(ur, 0.0, wfl.weight, wfr.weight);
				Vector2D p1 = ul.intersectAtTimeWeighted(ur, 1.0, wfl.weight, wfr.weight);
				if (p0 == null || p1 == null) { bi = new Vector2D(0, 0); pos0 = skNode.pos; }
				else { pos0 = p0; bi = p1.subtract(p0); }
			}
		}

		kv.velocity = bi;
		if (bi.getX() == 0 && bi.getY() == 0) { kv.infFast = true; kv.origin = skNode.pos; } else {
			kv.origin = pos0;
		}
		kv.ul = ul; kv.ur = ur; kv.wfl = wfl; kv.wfr = wfr;
		return kv;
	}

	public static void replaceInQueue(KineticTriangle t, double now, EventQueue queue, Deque<Event> imm) {
		if (t.event != null) { queue.discard(t.event); imm.remove(t.event); }
		Event e = CollapseEventComputer.compute(t, now, false);
		if (e != null) {
			queue.add(e);
		}
	}

	public static List<KineticTriangle> replaceKVertex(KineticTriangle t, KineticVertex v, KineticVertex newv, double now, IntUnaryOperator dir, EventQueue queue, Deque<Event> imm) {
		List<KineticTriangle> fan = new ArrayList<>();
		while (t != null) {
			int side = t.indexOfVertex(v);
			fan.add(t);
			t.vertices[side] = newv;
			if (newv.infFast && t.event != null) { queue.discard(t.event); imm.remove(t.event); } else {
				replaceInQueue(t, now, queue, imm);
			}
			t = t.neighbours[dir.applyAsInt(side)];
		}
		return fan;
	}

	public static void schedImm(KineticTriangle tri, double now, EventQueue queue, Deque<Event> imm) {
		if (tri.event != null) { queue.discard(tri.event); imm.remove(tri.event); }
		Event e = CollapseEventComputer.computeNewEdgeCollapse(tri, now);
		if (tri.getType() == 3) {
			e = new Event(now, tri, List.of(0,1,2), EventType.EDGE, 3);
		}
		tri.event = e;
		imm.addLast(e);
	}

	/** Equivalent of python update_circ(): update left/right links if kinetic. */
	public static void updateCirc(VertexRef left, VertexRef right, double now) {
		if (left instanceof KineticVertex) {
			((KineticVertex) left).setRight(right, now);
		}
		if (right instanceof KineticVertex) {
			((KineticVertex) right).setLeft(left, now);
		}
	}

	public static void flip(KineticTriangle t0, int side0, KineticTriangle t1, int side1) {
		int ax0 = apex(side0), or0 = orig(side0), de0 = dest(side0);
		int ax1 = apex(side1), or1 = orig(side1), de1 = dest(side1);
		VertexRef A = t0.vertices[ax0], B = t0.vertices[or0], C = t1.vertices[ax1], D = t0.vertices[de0];
		KineticTriangle AB = t0.neighbours[de0], BC = t1.neighbours[or1], CD = t1.neighbours[de1], DA = t0.neighbours[or0];

		KineticTriangle[] nghs = {AB, BC, CD, DA};
		VertexRef[] cors = {A, B, C, D};
		int[] sides = new int[4];
		for (int i=0; i<4; i++) {
			sides[i] = nghs[i] == null ? -1 : ccw(nghs[i].indexOfVertex(cors[i]));
		}

		KineticTriangle[] ts = {t0, t0, t1, t1};
		for (int i=0; i<4; i++) {
			if (nghs[i] != null) {
				nghs[i].neighbours[sides[i]] = ts[i];
			}
		}

		t0.vertices[0]=A; t0.vertices[1]=B; t0.vertices[2]=C; t0.neighbours[0]=BC; t0.neighbours[1]=t1; t0.neighbours[2]=AB;
		t1.vertices[0]=C; t1.vertices[1]=D; t1.vertices[2]=A; t1.neighbours[0]=DA; t1.neighbours[1]=t0; t1.neighbours[2]=CD;
	}

	/** Python handle_edge_event_1side: for type-3 triangle where only 1 side is reported collapsing. */
	public static void handleEdge1Side(Event evt, int step, Skeleton skel, EventQueue q, Deque<Event> imm) {
		KineticTriangle t = evt.tri;
		int e = evt.side.get(0);
		double now = evt.time;

		KineticVertex v0 = (KineticVertex) t.vertices[e];
		KineticVertex v1 = (KineticVertex) t.vertices[ccw(e)];
		KineticVertex v2 = (KineticVertex) t.vertices[cw(e)];

		stopKVertices(List.of(v1, v2), step, now, skel, null);
		KineticVertex kv = computeNewKVertex(v1.wfl, v2.wfr, now, v1.stopNode, skel.vertices.size()+1, v1.internal || v2.internal);
		skel.vertices.add(kv);

		stopKVertices(List.of(v0, kv), step, now, skel, null);
		t.stopsAt = now;
	}

	public static void handleEdge(Event evt, int step, Skeleton skel, EventQueue q, Deque<Event> imm) {
		KineticTriangle t = evt.tri;
		int e = evt.side.get(0);
		double now = evt.time;
		KineticVertex v1 = (KineticVertex)t.vertices[ccw(e)], v2 = (KineticVertex)t.vertices[cw(e)];

		// Use intersection of wavefront support lines at time now when possible.
		Vector2D pos = Geom.intersectionAtTimeWeighted(v1.wfl, v2.wfr, now);

		stopKVertices(List.of(v1, v2), step, now, skel, pos);
		KineticVertex kv = computeNewKVertex(v1.wfl, v2.wfr, now, v1.stopNode, skel.vertices.size()+1, v1.internal || v2.internal);
		skel.vertices.add(kv);

		VertexRef lRef = v1.getLeft(), rRef = v2.getRight();
		updateCirc(lRef, kv, now);
		updateCirc(kv, rRef, now);

		KineticTriangle aTri = t.neighbours[ccw(e)], bTri = t.neighbours[cw(e)], nTri = t.neighbours[e];
		List<KineticTriangle> fanA = List.of(), fanB = List.of();

		if (aTri != null) {
			aTri.neighbours[aTri.indexOfNeighbour(t)] = bTri;
			fanA = replaceKVertex(aTri, v2, kv, now, Core::cw, q, imm);
			if (!fanA.isEmpty()) {
				KineticTriangle tLast = fanA.get(fanA.size()-1);
				int sideIdx = cw(tLast.indexOfVertex(kv));
				VertexRef orig = tLast.vertices[ccw(sideIdx)];
				VertexRef dest = tLast.vertices[cw(sideIdx)];
				if (nearZero(Math.sqrt(orig.distance2At(dest, now)))) {
					schedImm(tLast, now, q, imm);
				}
			}
		}
		if (bTri != null) {
			bTri.neighbours[bTri.indexOfNeighbour(t)] = aTri;
			fanB = replaceKVertex(bTri, v1, kv, now, Core::ccw, q, imm);
			if (!fanB.isEmpty()) {
				KineticTriangle tLast = fanB.get(fanB.size()-1);
				int sideIdx = cw(tLast.indexOfVertex(kv));
				VertexRef orig = tLast.vertices[ccw(sideIdx)];
				VertexRef dest = tLast.vertices[cw(sideIdx)];
				if (nearZero(Math.sqrt(orig.distance2At(dest, now)))) {
					schedImm(tLast, now, q, imm);
				}
			}
		}
		if (nTri != null) {
			nTri.neighbours[nTri.indexOfNeighbour(t)] = null;
			if (nTri.event != null && nTri.stopsAt == null) {
				schedImm(nTri, now, q, imm);
			}
		}
		t.stopsAt = now;

		// Parallel/infFast special handling.
		if (kv.infFast) {
			if (!fanA.isEmpty() && !fanB.isEmpty()) {
				ArrayList<KineticTriangle> fan = new ArrayList<>(fanA);
				Collections.reverse(fan);
				fan.addAll(fanB);
				handleParallelFan(fan, kv, now, Core::ccw, step, skel, q, imm);
			} else if (!fanA.isEmpty()) {
				handleParallelFan(new ArrayList<>(fanA), kv, now, Core::cw, step, skel, q, imm);
			} else if (!fanB.isEmpty()) {
				handleParallelFan(new ArrayList<>(fanB), kv, now, Core::ccw, step, skel, q, imm);
			}
		}
	}

	public static void handleFlip(Event evt, double now, EventQueue q, Deque<Event> imm) {
		KineticTriangle t = evt.tri;
		int tSide = evt.side.get(0);
		KineticTriangle n = t.neighbours[tSide];
		int nSide = n.indexOfNeighbour(t);
		flip(t, tSide, n, nSide);
		replaceInQueue(t, now, q, imm);
		replaceInQueue(n, now, q, imm);
	}

	public static void handleSplit(Event evt, int step, Skeleton skel, EventQueue q, Deque<Event> imm) {
		KineticTriangle t = evt.tri;
		int e = evt.side.get(0);
		double now = evt.time;

		if (t.neighbours[e] != null) {
			return;
		}

		KineticVertex v = (KineticVertex) t.vertices[e];
		KineticVertex v1 = (KineticVertex) t.vertices[ccw(e)];
		KineticVertex v2 = (KineticVertex) t.vertices[cw(e)];

		stopKVertices(List.of(v), step, now, skel, null);
		SkeletonNode skNode = v.stopNode;

		KineticVertex vb = computeNewKVertex(v.wfl, v2.wfl, now, skNode, skel.vertices.size()+1, v.internal || v2.internal);
		skel.vertices.add(vb);
		KineticVertex va = computeNewKVertex(v1.wfr, v.wfr, now, skNode, skel.vertices.size()+1, v.internal || v1.internal);
		skel.vertices.add(va);

		updateCirc(v.getLeft(), vb, now);
		updateCirc(vb, v2, now);
		updateCirc(v1, va, now);
		updateCirc(va, v.getRight(), now);

		KineticTriangle b = t.neighbours[ccw(e)];
		if (b != null) {
			b.neighbours[b.indexOfNeighbour(t)] = null;
		}
		List<KineticTriangle> fanB = replaceKVertex(b, v, vb, now, Core::ccw, q, imm);

		KineticTriangle a = t.neighbours[cw(e)];
		if (a != null) {
			a.neighbours[a.indexOfNeighbour(t)] = null;
		}
		List<KineticTriangle> fanA = replaceKVertex(a, v, va, now, Core::cw, q, imm);

		t.stopsAt = now;

		if (va.infFast && !fanA.isEmpty()) {
			handleParallelFan(new ArrayList<>(fanA), va, now, Core::cw, step, skel, q, imm);
		}
		if (vb.infFast && !fanB.isEmpty()) {
			handleParallelFan(new ArrayList<>(fanB), vb, now, Core::ccw, step, skel, q, imm);
		}
	}

	private static void handleParallelFan(
			List<KineticTriangle> fan,
			KineticVertex pivot,
			double now,
			IntUnaryOperator direction,
			int step,
			Skeleton skel,
			EventQueue q,
			Deque<Event> imm
			) {
		if (fan == null || fan.isEmpty()) {
			throw new IllegalArgumentException("expected fan of triangles");
		}
		if (!pivot.infFast) {
			return;
		}

		KineticTriangle firstTri = fan.get(0);
		if (firstTri.getType() == 3) {
			double[] dists = new double[3];
			for (int side=0; side<3; side++) {
				VertexRef s = firstTri.vertices[ccw(side)];
				VertexRef e = firstTri.vertices[cw(side)];
				dists[side] = s.positionAt(now).distance(e.positionAt(now));
			}
			double mn = Math.min(dists[0], Math.min(dists[1], dists[2]));
			int ct = 0, idx = -1;
			for (int i=0; i<3; i++) {
				if (nearZero(dists[i] - mn)) { ct++; idx = i; }
			}
			if (nearZero(mn) && ct == 1) {
				KineticVertex pivot2 = (KineticVertex) firstTri.vertices[idx];
				handleParallelEdgeEventEvenLegs(firstTri, idx, pivot2, now, step, skel, q, imm);
				return;
			}
			handleParallelEdgeEvent3Tri(firstTri, firstTri.indexOfVertex(pivot), pivot, now, step, skel, q, imm);
			return;
		}

		int d0 = direction.applyAsInt(0);
		final boolean isCw;
		if (d0 == 2) {
			isCw = true;
		} else if (d0 == 1) {
			isCw = false;
		} else {
			throw new IllegalStateException("Unexpected direction operator: op(0)=" + d0);
		}
		KineticTriangle left = isCw ? fan.get(0) : fan.get(fan.size()-1);
		KineticTriangle right = isCw ? fan.get(fan.size()-1) : fan.get(0);

		int leftLegIdx = ccw(left.indexOfVertex(pivot));
		VertexRef vls = left.vertices[ccw(leftLegIdx)];
		VertexRef vle = left.vertices[cw(leftLegIdx)];
		double leftDist = vls.positionAt(now).distance(vle.positionAt(now));

		int rightLegIdx = cw(right.indexOfVertex(pivot));
		VertexRef vrs = right.vertices[ccw(rightLegIdx)];
		VertexRef vre = right.vertices[cw(rightLegIdx)];
		double rightDist = vrs.positionAt(now).distance(vre.positionAt(now));

		double mn = Math.min(leftDist, rightDist);
		boolean leftMin = nearZero(leftDist - mn);
		boolean rightMin = nearZero(rightDist - mn);

		if (leftMin && rightMin) {
			if (fan.size() == 1) {
				handleParallelEdgeEventEvenLegs(firstTri, firstTri.indexOfVertex(pivot), pivot, now, step, skel, q, imm);
			} else if (fan.size() == 2) {
				boolean all2 = true;
				for (KineticTriangle t : fan) {
					int lIdx = ccw(t.indexOfVertex(pivot));
					int rIdx = cw(t.indexOfVertex(pivot));
					VertexRef l1 = t.vertices[ccw(lIdx)], l2 = t.vertices[cw(lIdx)];
					VertexRef r1 = t.vertices[ccw(rIdx)], r2 = t.vertices[cw(rIdx)];
					double ld = l1.positionAt(now).distance(l2.positionAt(now));
					double rd = r1.positionAt(now).distance(r2.positionAt(now));
					double mm = Math.min(ld, rd);
					int u = 0;
					if (nearZero(ld - mm)) {
						u++;
					}
					if (nearZero(rd - mm)) {
						u++;
					}
					if (u != 2) {
						all2 = false;
					}
				}
				if (all2) {
					for (KineticTriangle t : fan) {
						handleParallelEdgeEventEvenLegs(t, t.indexOfVertex(pivot), pivot, now, step, skel, q, imm);
					}
				} else {
					KineticTriangle t0 = fan.get(0), t1 = fan.get(1);
					int side0 = t0.indexOfNeighbour(t1);
					int side1 = t1.indexOfNeighbour(t0);
					flip(t0, side0, t1, side1);
					if (hasInfFast(t0)) {
						handleParallelEdgeEventEvenLegs(t0, t0.indexOfVertex(pivot), pivot, now, step, skel, q, imm);
					}
					if (hasInfFast(t1)) {
						handleParallelEdgeEventEvenLegs(t1, t1.indexOfVertex(pivot), pivot, now, step, skel, q, imm);
					}
				}
			} else {
				throw new UnsupportedOperationException("More than 2 triangles in equal-legs parallel fan");
			}
			return;
		}

		if (rightMin) {
			handleParallelEdgeEventShorterLeg(right, rightLegIdx, pivot, now, step, skel, q, imm);
		} else {
			handleParallelEdgeEventShorterLeg(left, leftLegIdx, pivot, now, step, skel, q, imm);
		}
	}

	private static boolean hasInfFast(KineticTriangle t) {
		for (VertexRef v : t.vertices) {
			if (v instanceof KineticVertex && ((KineticVertex) v).infFast) {
				return true;
			}
		}
		return false;
	}

	private static void handleParallelEdgeEventShorterLeg(
			KineticTriangle t, int e, KineticVertex pivot, double now, int step,
			Skeleton skel, EventQueue q, Deque<Event> imm
			) {
		KineticVertex v1 = (KineticVertex) t.vertices[ccw(e)];
		KineticVertex v2 = (KineticVertex) t.vertices[cw(e)];

		List<KineticVertex> toStop = new ArrayList<>();
		if (!v1.infFast) {
			toStop.add(v1);
		}
		if (!v2.infFast) {
			toStop.add(v2);
		}
		if (!toStop.isEmpty()) {
			stopKVertices(toStop, step, now, skel, null);
		}

		SkeletonNode skNode = (!toStop.isEmpty()) ? toStop.get(0).stopNode : (pivot.stopNode != null ? pivot.stopNode : pivot.startNode);
		if (skNode == null) {
			return;
		}

		if (pivot.stopNode == null) { pivot.stopNode = skNode; pivot.stopsAt = now; }
		t.stopsAt = now;

		KineticVertex kv = computeNewKVertex(v1.wfl, v2.wfr, now, skNode, skel.vertices.size()+1, v1.internal || v2.internal);
		// Python: overwrite wfl/wfr from the circular neighbors (or None).
		VertexRef leftRef = v1.getLeft();
		if (leftRef == null) {
			kv.wfl = null;
		} else if (leftRef instanceof KineticVertex) {
			kv.wfl = ((KineticVertex) leftRef).wfr;
		} else {
			throw new IllegalStateException("Expected v1.left to be KineticVertex or null");
		}

		VertexRef rightRef = v2.getRight();
		if (rightRef == null) {
			kv.wfr = null;
		} else if (rightRef instanceof KineticVertex) {
			kv.wfr = ((KineticVertex) rightRef).wfl;
		} else {
			throw new IllegalStateException("Expected v2.right to be KineticVertex or null");
		}
		skel.vertices.add(kv);

		updateCirc(v1.getLeft(), kv, now);
		updateCirc(kv, v2.getRight(), now);

		KineticTriangle a = t.neighbours[ccw(e)];
		KineticTriangle b = t.neighbours[cw(e)];
		KineticTriangle n = t.neighbours[e];

		List<KineticTriangle> fanA = List.of(), fanB = List.of();
		if (a != null) {
			a.neighbours[a.indexOfNeighbour(t)] = b;
			fanA = replaceKVertex(a, v2, kv, now, Core::cw, q, imm);
		}
		if (b != null) {
			b.neighbours[b.indexOfNeighbour(t)] = a;
			fanB = replaceKVertex(b, v1, kv, now, Core::ccw, q, imm);
		}
		if (n != null) {
			n.neighbours[n.indexOfNeighbour(t)] = null;
			if (n.event != null && n.stopsAt == null) {
				schedImm(n, now, q, imm);
			}
		}

		if (kv.infFast) {
			if (!fanA.isEmpty() && !fanB.isEmpty()) {
				ArrayList<KineticTriangle> fan = new ArrayList<>(fanA);
				Collections.reverse(fan);
				fan.addAll(fanB);
				handleParallelFan(fan, kv, now, Core::ccw, step, skel, q, imm);
			} else if (!fanA.isEmpty()) {
				handleParallelFan(new ArrayList<>(fanA), kv, now, Core::cw, step, skel, q, imm);
			} else if (!fanB.isEmpty()) {
				handleParallelFan(new ArrayList<>(fanB), kv, now, Core::ccw, step, skel, q, imm);
			}
		}
	}

	private static void handleParallelEdgeEventEvenLegs(
			KineticTriangle t, int e, KineticVertex pivot, double now, int step,
			Skeleton skel, EventQueue q, Deque<Event> imm
			) {
		KineticVertex v1 = (KineticVertex) t.vertices[ccw(e)];
		KineticVertex v2 = (KineticVertex) t.vertices[cw(e)];
		stopKVertices(List.of(v1, v2), step, now, skel, null);
		pivot.stopNode = v1.stopNode;
		pivot.stopsAt = now;
		t.stopsAt = now;

		KineticTriangle n = t.neighbours[e];
		if (n != null) {
			n.neighbours[n.indexOfNeighbour(t)] = null;
			if (n.event != null && n.stopsAt == null) {
				schedImm(n, now, q, imm);
			}
		}
	}

	private static void handleParallelEdgeEvent3Tri(
			KineticTriangle t, int e, KineticVertex pivot, double now, int step,
			Skeleton skel, EventQueue q, Deque<Event> imm
			) {
		KineticVertex v1 = (KineticVertex) t.vertices[ccw(e)];
		KineticVertex v2 = (KineticVertex) t.vertices[cw(e)];

		double m1 = v1.velocityAt(now).length();
		double m2 = v2.velocityAt(now).length();

		if (m2 < m1) {
			stopKVertices(List.of(v2), step, now, skel, null);
			v1.stopNode = v2.stopNode;
			v1.stopsAt = now;
		} else {
			stopKVertices(List.of(v1), step, now, skel, null);
			v2.stopNode = v1.stopNode;
			v2.stopsAt = now;
		}

		pivot.stopNode = v1.stopNode != null ? v1.stopNode : v2.stopNode;
		pivot.stopsAt = now;
		t.stopsAt = now;
	}

	public static double eventLoop(EventQueue q, Skeleton skel) {
		double now = 0;
		int step = 0;
		Deque<Event> imm = new ArrayDeque<>();
		int guard = 0;

		while (!q.isEmpty() || !imm.isEmpty()) {
			if (guard++ > LOOP_MAX) {
				throw new RuntimeException("Exceeded maximum event loop iterations (" + LOOP_MAX + ") (probably degenerate input in some way that is unhandled)");
			}
			step++;
			Event evt = imm.isEmpty() ? q.pop() : imm.pollFirst();
			if (evt == null) {
				continue;
			}
			now = evt.time;

			if (evt.tri.stopsAt != null) {
				continue;
			}

			if (evt.tp == EventType.EDGE) {
				if (evt.side.size() == 3) {
					stopKVertices(Arrays.asList((KineticVertex)evt.tri.vertices[0], (KineticVertex)evt.tri.vertices[1], (KineticVertex)evt.tri.vertices[2]), step, now, skel, null);
					for (var n : evt.tri.neighbours) {
						if (n != null && n.event != null && n.stopsAt == null) {
							n.neighbours[n.indexOfNeighbour(evt.tri)] = null;
							schedImm(n, now, q, imm);
						}
					}
					evt.tri.stopsAt = now;
				} else if (evt.side.size() == 2) {
					throw new RuntimeException("Impossible configuration: triangle has 2 sides collapsing");
				} else if (evt.side.size() == 1 && evt.tri.getType() == 3) {
					handleEdge1Side(evt, step, skel, q, imm);
				} else {
					handleEdge(evt, step, skel, q, imm);
				}
			} else if (evt.tp == EventType.FLIP) {
				handleFlip(evt, now, q, imm);
			} else if (evt.tp == EventType.SPLIT) {
				handleSplit(evt, step, skel, q, imm);
			} else {
				throw new RuntimeException("Unknown event type: " + evt.tp);
			}
		}
		return now;
	}
}
