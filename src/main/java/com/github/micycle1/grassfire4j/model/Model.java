package com.github.micycle1.grassfire4j.model;

import static com.github.micycle1.grassfire4j.geom.Geom.STOP_EPS;
import static com.github.micycle1.grassfire4j.geom.Geom.dist2;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.math.Vector2D;
import org.locationtech.jts.operation.polygonize.Polygonizer;

import com.github.micycle1.grassfire4j.geom.Geom.Line2;
import com.github.micycle1.grassfire4j.geom.Geom.WaveFront;

/**
 * Domain model types for the kinetic straight-skeleton computation.
 * <p>
 * These types represent solver state after input has been normalised into
 * {@code InputMesh}, independent of any specific adapter implementation.
 */
public final class Model {

	public static class SkeletonNode {
		public final Vector2D pos;
		public final int step;
		public final Integer info;
		public SkeletonNode(Vector2D pos, int step, Integer info) { this.pos = pos; this.step = step; this.info = info; }
	}

	public interface VertexRef {
		Vector2D positionAt(double time);
		Vector2D velocityAt(double time);
		default double distance2At(VertexRef other, double time) {
			return dist2(positionAt(time), other.positionAt(time));
		}
	}

	public static class InfiniteVertex implements VertexRef {
		public Vector2D origin;
		public boolean internal = false;
		public int info = 0;
		public InfiniteVertex(Vector2D origin) { this.origin = origin; }
		@Override public Vector2D positionAt(double time) { return origin; }
		@Override public Vector2D velocityAt(double time) { return new Vector2D(0, 0); }
	}

	public static class KineticVertex implements VertexRef {
		public enum Turn {
			RIGHT_REFLEX,
			LEFT_CONVEX,
			STRAIGHT
		}

		public Vector2D origin, velocity;
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

		@Override public Vector2D positionAt(double time) {
			if (infFast) {
				return startNode.pos;
			}
			if (stopsAt != null && stopNode != null && time >= stopsAt - STOP_EPS) {
				return stopNode.pos;
			}
			return new Vector2D(origin.getX() + time * velocity.getX(), origin.getY() + time * velocity.getY());
		}

		@Override public Vector2D velocityAt(double time) {
			if (infFast || (stopsAt != null && time >= stopsAt - STOP_EPS)) {
				return new Vector2D(0, 0);
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
		public final int sideMask;
		public final EventType tp;
		public final int triangleTp;
		public boolean valid = true;
		public long counter;

		public Event(double time, KineticTriangle tri, int sideMask, EventType tp, int triangleTp) {
			this.time = time; this.tri = tri; this.sideMask = sideMask; this.tp = tp; this.triangleTp = triangleTp;
		}

		public int sideCount() { return Integer.bitCount(sideMask); }

		public int singleSide() {
			if (sideCount() != 1) {
				throw new IllegalStateException("Event does not have exactly one side (mask=" + sideMask + ")");
			}
			return Integer.numberOfTrailingZeros(sideMask);
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

	/**
	 * Straight-skeleton result model produced by the kinetic solver.
	 * <p>
	 * Contains skeleton nodes, kinetic vertices, and helper methods to export the
	 * resulting linework and derived polygon faces in JTS geometry form.
	 */
	public static class Skeleton {
		public List<SkeletonNode> skNodes = new ArrayList<>();
		public List<KineticVertex> vertices = new ArrayList<>();
		public List<KineticTriangle> triangles = new ArrayList<>();

		public record Segment(Coordinate p1, Coordinate p2, Integer info1, Integer info2) {}

		/**
		 * Returns all skeleton segments as endpoint pairs.
		 * <p>
		 * Stopped vertices produce finite segments from start node to stop node.
		 * Non-stopped vertices are represented by a long ray sampled at {@code t=1000}.
		 *
		 * @return list of skeleton segments
		 */
		public List<Segment> segments() {
			List<Segment> segs = new ArrayList<>();
			for (var v : vertices) {
				if (v.stopsAt != null) {
					if (v.startNode == v.stopNode) {
						continue;
					}
					segs.add(new Segment(v.startNode.pos.toCoordinate(), v.stopNode.pos.toCoordinate(), v.startNode.info, v.stopNode.info));
				} else {
					segs.add(new Segment(v.startNode.pos.toCoordinate(), v.positionAt(1000.0).toCoordinate(), v.startNode.info, null));
				}
			}
			return segs;
		}

		/**
		 * Returns skeleton segments as a JTS {@link MultiLineString} using a default
		 * {@link GeometryFactory}.
		 *
		 * @return skeleton linework
		 */
		public MultiLineString asMultiLineString() {
			return asMultiLineString(new GeometryFactory());
		}

		/**
		 * Returns skeleton segments as a JTS {@link MultiLineString}.
		 *
		 * @param geometryFactory geometry factory used to build the output geometry
		 * @return skeleton linework
		 */
		public MultiLineString asMultiLineString(GeometryFactory geometryFactory) {
			List<Segment> segs = segments();
			LineString[] lines = new LineString[segs.size()];
			for (int i = 0; i < segs.size(); i++) {
				Segment s = segs.get(i);
				lines[i] = geometryFactory.createLineString(new Coordinate[] { s.p1(), s.p2() });
			}
			return geometryFactory.createMultiLineString(lines);
		}

		/**
		 * Polygonizes the partition induced by the input polygon boundary and skeleton
		 * linework.
		 * <p>
		 * The caller should pass the same polygon used to compute this skeleton.
		 *
		 * @param polygon source polygon used for skeleton construction
		 * @return polygonized face geometry
		 */
		public Geometry asPolygonFaces(Polygon polygon) {
			Objects.requireNonNull(polygon, "polygon");
			Geometry coverage = polygon.getBoundary().union(asMultiLineString(polygon.getFactory()));
			Polygonizer polygonizer = new Polygonizer();
			polygonizer.setCheckRingsValid(false);
			polygonizer.add(coverage);
			return polygonizer.getGeometry();
		}
	}
}
