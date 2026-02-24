package com.github.micycle1.grassfire4j.model;

import static com.github.micycle1.grassfire4j.geom.Geom.STOP_EPS;
import static com.github.micycle1.grassfire4j.geom.Geom.dist2;

import java.util.ArrayList;
import java.util.List;

import com.github.micycle1.grassfire4j.geom.Geom.Line2;
import com.github.micycle1.grassfire4j.geom.Geom.Vec2;
import com.github.micycle1.grassfire4j.geom.Geom.WaveFront;

/**
 * Domain model types for the kinetic straight-skeleton computation.
 * <p>
 * These types represent solver state after input has been normalised into
 * {@code InputMesh}, independent of any specific adapter implementation.
 */
public final class Model {

	public static class SkeletonNode {
		public final Vec2 pos;
		public final int step;
		public final Integer info;
		public SkeletonNode(Vec2 pos, int step, Integer info) { this.pos = pos; this.step = step; this.info = info; }
	}

	public interface VertexRef {
		Vec2 positionAt(double time);
		Vec2 velocityAt(double time);
		default double distance2At(VertexRef other, double time) {
			return dist2(positionAt(time), other.positionAt(time));
		}
	}

	public static class InfiniteVertex implements VertexRef {
		public Vec2 origin;
		public boolean internal = false;
		public int info = 0;
		public InfiniteVertex(Vec2 origin) { this.origin = origin; }
		@Override public Vec2 positionAt(double time) { return origin; }
		@Override public Vec2 velocityAt(double time) { return new Vec2(0, 0); }
	}

	public static class KineticVertex implements VertexRef {
		public enum Turn {
			RIGHT_REFLEX,
			LEFT_CONVEX,
			STRAIGHT
		}

		public Vec2 origin, velocity;
		public Double startsAt, stopsAt;
		public SkeletonNode startNode, stopNode;
		public Line2 ul, ur;
		public WaveFront wfl, wfr;
		public int info = 0;
		public boolean infFast = false, internal = false;
		public Turn turn;

		public record HistEntry(double start, Double stop, VertexRef ref) {}
		public final List<HistEntry> leftHist = new ArrayList<>();
		public final List<HistEntry> rightHist = new ArrayList<>();

		public boolean isStopped() { return stopNode != null; }

		@Override public Vec2 positionAt(double time) {
			if (infFast) {
				return startNode.pos;
			}
			if (stopsAt != null && stopNode != null && time >= stopsAt - STOP_EPS) {
				return stopNode.pos;
			}
			return new Vec2(origin.x() + time * velocity.x(), origin.y() + time * velocity.y());
		}

		@Override public Vec2 velocityAt(double time) {
			if (infFast || (stopsAt != null && time >= stopsAt - STOP_EPS)) {
				return new Vec2(0, 0);
			}
			return velocity;
		}

		public VertexRef getLeft() { return leftHist.isEmpty() ? null : leftHist.get(leftHist.size()-1).ref(); }
		public VertexRef getRight() { return rightHist.isEmpty() ? null : rightHist.get(rightHist.size()-1).ref(); }

		public void setLeft(VertexRef ref, double time) {
			if (!leftHist.isEmpty()) {
				var last = leftHist.get(leftHist.size()-1);
				leftHist.set(leftHist.size()-1, new HistEntry(last.start(), time, last.ref()));
			}
			leftHist.add(new HistEntry(time, null, ref));
		}

		public void setRight(VertexRef ref, double time) {
			if (!rightHist.isEmpty()) {
				var last = rightHist.get(rightHist.size()-1);
				rightHist.set(rightHist.size()-1, new HistEntry(last.start(), time, last.ref()));
			}
			rightHist.add(new HistEntry(time, null, ref));
		}
	}

	public static class KineticTriangle {
		public final VertexRef[] vertices = new VertexRef[3];
		public final KineticTriangle[] neighbours = new KineticTriangle[3];
		public final WaveFront[] wavefrontSupportLines = new WaveFront[3];
		public Event event;
		public int info = 0, uid = 0;
		public Double stopsAt;
		public boolean internal = false;

		public int getType() {
			int noneCt = 0;
			for (var n : neighbours) {
				if (n == null) {
					noneCt++;
				}
			}
			return noneCt;
		}

		public boolean isFinite() {
			for (var v : vertices) {
				if (!(v instanceof KineticVertex)) {
					return false;
				}
			}
			return true;
		}

		public int indexOfVertex(VertexRef v) {
			for (int i = 0; i < 3; i++) {
				if (vertices[i] == v) {
					return i;
				}
			}
			return -1;
		}

		public int indexOfNeighbour(KineticTriangle t) {
			for (int i = 0; i < 3; i++) {
				if (neighbours[i] == t) {
					return i;
				}
			}
			return -1;
		}
	}

	public static class Event implements Comparable<Event> {
		public enum EventType {
			EDGE,
			FLIP,
			SPLIT
		}

		public final double time;
		public final KineticTriangle tri;
		public final List<Integer> side;
		public final EventType tp;
		public final int triangleTp;
		public boolean valid = true;
		public long counter;

		public Event(double time, KineticTriangle tri, List<Integer> side, EventType tp, int triangleTp) {
			this.time = time; this.tri = tri; this.side = side; this.tp = tp; this.triangleTp = triangleTp;
		}

		@Override public int compareTo(Event o) {
			int cmp = Double.compare(this.time, o.time);
			if (cmp != 0) {
				return cmp;
			}
			cmp = Integer.compare(o.triangleTp, this.triangleTp); // Descending
			if (cmp != 0) {
				return cmp;
			}
			cmp = Integer.compare(this.tri.uid, o.tri.uid);
			if (cmp != 0) {
				return cmp;
			}
			return Long.compare(this.counter, o.counter);
		}
	}

	public static class Skeleton {
		public List<SkeletonNode> skNodes = new ArrayList<>();
		public List<KineticVertex> vertices = new ArrayList<>();
		public List<KineticTriangle> triangles = new ArrayList<>();

		public record Segment(Vec2 p1, Vec2 p2, Integer info1, Integer info2) {}
		public List<Segment> segments() {
			List<Segment> segs = new ArrayList<>();
			for (var v : vertices) {
				if (v.stopsAt != null) {
					if (v.startNode == v.stopNode) {
						continue;
					}
					segs.add(new Segment(v.startNode.pos, v.stopNode.pos, v.startNode.info, v.stopNode.info));
				} else {
					segs.add(new Segment(v.startNode.pos, v.positionAt(1000.0), v.startNode.info, null));
				}
			}
			return segs;
		}
	}
}
