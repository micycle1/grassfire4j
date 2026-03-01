package com.github.micycle1.grassfire4j.model;

import static com.github.micycle1.grassfire4j.geom.Geom.STOP_EPS;
import static com.github.micycle1.grassfire4j.geom.Geom.dist2;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
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

	private Model() {
	}

	public static class SkeletonNode {
		public final Vector2D pos;
		public final int step;
		public final Integer info;

		public SkeletonNode(Vector2D pos, int step, Integer info) {
			this.pos = pos;
			this.step = step;
			this.info = info;
		}
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

		public InfiniteVertex(Vector2D origin) {
			this.origin = origin;
		}

		@Override
		public Vector2D positionAt(double time) {
			return origin;
		}

		@Override
		public Vector2D velocityAt(double time) {
			return new Vector2D(0, 0);
		}
	}

	public static class KineticVertex implements VertexRef {
		public enum Turn {
			RIGHT_REFLEX, LEFT_CONVEX, STRAIGHT
		}

		public Vector2D origin, velocity;
		public Double startsAt, stopsAt;
		public SkeletonNode startNode, stopNode;
		public Line2 ul, ur;
		public WaveFront wfl, wfr;
		public int info = 0;
		public boolean infFast = false, internal = false;
		public Turn turn;

		public record HistEntry(double start, Double stop, VertexRef ref) {
		}

		public final List<HistEntry> leftHist = new ArrayList<>();
		public final List<HistEntry> rightHist = new ArrayList<>();

		public boolean isStopped() {
			return stopNode != null;
		}

		@Override
		public Vector2D positionAt(double time) {
			if (infFast) {
				return startNode.pos;
			}
			if (stopsAt != null && stopNode != null && time >= stopsAt - STOP_EPS) {
				return stopNode.pos;
			}
			return new Vector2D(origin.getX() + time * velocity.getX(), origin.getY() + time * velocity.getY());
		}

		@Override
		public Vector2D velocityAt(double time) {
			if (infFast || (stopsAt != null && time >= stopsAt - STOP_EPS)) {
				return new Vector2D(0, 0);
			}
			return velocity;
		}

		public VertexRef getLeft() {
			return leftHist.isEmpty() ? null : leftHist.get(leftHist.size() - 1).ref();
		}

		public VertexRef getRight() {
			return rightHist.isEmpty() ? null : rightHist.get(rightHist.size() - 1).ref();
		}

		public void setLeft(VertexRef ref, double time) {
			if (!leftHist.isEmpty()) {
				var last = leftHist.get(leftHist.size() - 1);
				leftHist.set(leftHist.size() - 1, new HistEntry(last.start(), time, last.ref()));
			}
			leftHist.add(new HistEntry(time, null, ref));
		}

		public void setRight(VertexRef ref, double time) {
			if (!rightHist.isEmpty()) {
				var last = rightHist.get(rightHist.size() - 1);
				rightHist.set(rightHist.size() - 1, new HistEntry(last.start(), time, last.ref()));
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
			EDGE, FLIP, SPLIT
		}

		public final double time;
		public final KineticTriangle tri;
		public final int sideMask;
		public final EventType tp;
		public final int triangleTp;
		public boolean valid = true;
		public long counter;

		public Event(double time, KineticTriangle tri, int sideMask, EventType tp, int triangleTp) {
			this.time = time;
			this.tri = tri;
			this.sideMask = sideMask;
			this.tp = tp;
			this.triangleTp = triangleTp;
		}

		public int sideCount() {
			return Integer.bitCount(sideMask);
		}

		public int singleSide() {
			if (sideCount() != 1) {
				throw new IllegalStateException("Event does not have exactly one side (mask=" + sideMask + ")");
			}
			return Integer.numberOfTrailingZeros(sideMask);
		}

		@Override
		public int compareTo(Event o) {
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
		private static final int OUTSIDE_FACE = -1;
		private static final double FACE_RING_AREA_EPS = 1e-12;

		private static record CoordKey(double x, double y) {
		}

		private static final class HalfEdge {
			final int from;
			final int to;
			final int twin;
			final int leftFace;
			final double angle;
			int next = -1;

			HalfEdge(int from, int to, int twin, int leftFace, double angle) {
				this.from = from;
				this.to = to;
				this.twin = twin;
				this.leftFace = leftFace;
				this.angle = angle;
			}
		}

		public List<SkeletonNode> skNodes = new ArrayList<>();
		public List<KineticVertex> vertices = new ArrayList<>();
		public List<KineticTriangle> triangles = new ArrayList<>();

		public record BoundaryEdge(Coordinate from, Coordinate to, int edgeId) {
		}

		public List<BoundaryEdge> boundaryEdges = List.of();

		public record Segment(Coordinate p1, Coordinate p2, Integer info1, Integer info2) {
		}

		/**
		 * Returns all skeleton segments as endpoint pairs.
		 * <p>
		 * Stopped vertices produce finite segments from start node to stop node.
		 *
		 * @return list of skeleton segments
		 */
		public List<Segment> segments() {
			List<Segment> segs = new ArrayList<>();
			for (var v : vertices) {
				double startTime = v.startsAt == null ? 0.0 : v.startsAt.doubleValue();
				if (v.stopsAt != null) {
					if (v.startNode == v.stopNode) {
						continue;
					}
					segs.add(new Segment(new Coordinate(v.startNode.pos.getX(), v.startNode.pos.getY(), startTime),
							new Coordinate(v.stopNode.pos.getX(), v.stopNode.pos.getY(), v.stopsAt.doubleValue()), v.startNode.info, v.stopNode.info));
				} else {
					Vector2D rayEnd = v.positionAt(1000.0);
					segs.add(new Segment(new Coordinate(v.startNode.pos.getX(), v.startNode.pos.getY(), startTime),
							new Coordinate(rayEnd.getX(), rayEnd.getY(), 1000.0), v.startNode.info, null));
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
		 * The caller should pass the same polygon used to compute this skeleton. Faces
		 * are derived from a half-edge walk over boundary edges and skeleton arcs
		 * labeled by source boundary edge IDs. If boundary-edge labels are unavailable
		 * (e.g. custom adapters not providing IDs), this method falls back to
		 * polygonizer-based extraction.
		 *
		 * @param polygon source polygon used for skeleton construction
		 * @return polygonized face geometry
		 */
		public Geometry asPolygonFaces() {
			Geometry g = buildFacesFromTopology();
			return g == null ? new GeometryFactory().createGeometryCollection() : g;
		}

		private Geometry buildFacesFromTopology() {
			GeometryFactory factory = new GeometryFactory();

			List<HalfEdge> halfEdges = new ArrayList<>();
			List<List<Integer>> outgoing = new ArrayList<>();
			List<Coordinate> nodes = new ArrayList<>();
			Map<CoordKey, Integer> nodeLookup = new HashMap<>();

			// 1) Add stored boundary edges (interior is on LEFT of from->to)
			for (BoundaryEdge e : boundaryEdges) {
				addEdgePair(e.from(), e.to(), e.edgeId(), // leftFaceForward
						null, // leftFaceBackward = OUTSIDE
						halfEdges, outgoing, nodes, nodeLookup);
			}

			// 2) Add skeleton edges
			boolean hasLabeledSkeletonEdges = false;
			for (var v : vertices) {
				if (v.stopsAt == null || v.startNode == null || v.stopNode == null || v.startNode == v.stopNode) {
					continue;
				}
				Integer leftFace = wavefrontFaceId(v.wfl);
				Integer rightFace = wavefrontFaceId(v.wfr);
				if (leftFace == null && rightFace == null) {
					continue;
				}
				hasLabeledSkeletonEdges = true;
				double startTime = v.startsAt == null ? 0.0 : v.startsAt.doubleValue();
				addEdgePair(new Coordinate(v.startNode.pos.getX(), v.startNode.pos.getY(), startTime),
						new Coordinate(v.stopNode.pos.getX(), v.stopNode.pos.getY(), v.stopsAt.doubleValue()), leftFace, rightFace, halfEdges, outgoing, nodes,
						nodeLookup);
			}

			if (!hasLabeledSkeletonEdges || halfEdges.isEmpty()) {
				return null;
			}

			// 3) Sort outgoing edges by angle (same as before)
			int[] indexInOutgoing = new int[halfEdges.size()];
			for (List<Integer> edgesAtNode : outgoing) {
				edgesAtNode.sort(Comparator.comparingDouble(idx -> halfEdges.get(idx).angle));
				for (int i = 0; i < edgesAtNode.size(); i++) {
					indexInOutgoing[edgesAtNode.get(i)] = i;
				}
			}

			// 4) Build next pointers (same as before)
			for (HalfEdge halfEdge : halfEdges) {
				HalfEdge e = halfEdge;
				List<Integer> aroundTo = outgoing.get(e.to);
				if (aroundTo.isEmpty()) {
					return null;
				}
				int twinPos = indexInOutgoing[e.twin];
				int nextPos = (twinPos - 1 + aroundTo.size()) % aroundTo.size();
				e.next = aroundTo.get(nextPos);
			}

			// 5) Extract faces by walking half-edges (same as before)
			boolean[] visited = new boolean[halfEdges.size()];
			List<Polygon> faces = new ArrayList<>();

			for (int i = 0; i < halfEdges.size(); i++) {
				HalfEdge start = halfEdges.get(i);
				if (start.leftFace == OUTSIDE_FACE || visited[i]) {
					continue;
				}

				List<Coordinate> ring = new ArrayList<>();
				int current = i;
				int guard = 0;

				while (true) {
					if (guard++ > halfEdges.size() + 5) {
						return null;
					}
					HalfEdge edge = halfEdges.get(current);
					if (edge.leftFace != start.leftFace || visited[current]) {
						break;
					}
					visited[current] = true;
					ring.add(nodes.get(edge.from));
					current = edge.next;
					if (current == i || current < 0) {
						break;
					}
				}

				Coordinate[] coords = sanitizeRing(ring);
				if (coords == null) {
					continue;
				}

				Polygon face = factory.createPolygon(coords);
				if (face.getArea() <= FACE_RING_AREA_EPS) {
					continue;
				}
				faces.add(face);
			}

			if (faces.isEmpty()) {
				return null;
			}
			return factory.buildGeometry(faces);
		}

		private static Integer wavefrontFaceId(WaveFront wavefront) {
			if (wavefront == null || !(wavefront.data instanceof Integer face)) {
				return null;
			}
			return face.intValue() >= 0 ? face : null;
		}

		private static void addEdgePair(Coordinate fromCoord, Coordinate toCoord, Integer leftFaceForward, Integer leftFaceBackward, List<HalfEdge> halfEdges,
				List<List<Integer>> outgoing, List<Coordinate> nodes, Map<CoordKey, Integer> nodeLookup) {
			int from = nodeIndex(fromCoord, nodes, nodeLookup, outgoing);
			int to = nodeIndex(toCoord, nodes, nodeLookup, outgoing);
			if (from == to) {
				return;
			}

			Coordinate a = nodes.get(from);
			Coordinate b = nodes.get(to);
			double abx = b.x - a.x;
			double aby = b.y - a.y;
			double bax = -abx;
			double bay = -aby;
			double angleAB = Math.atan2(aby, abx);
			double angleBA = Math.atan2(bay, bax);
			int abIdx = halfEdges.size();
			int baIdx = abIdx + 1;
			halfEdges.add(new HalfEdge(from, to, baIdx, leftFaceForward == null ? OUTSIDE_FACE : leftFaceForward.intValue(), angleAB));
			halfEdges.add(new HalfEdge(to, from, abIdx, leftFaceBackward == null ? OUTSIDE_FACE : leftFaceBackward.intValue(), angleBA));
			outgoing.get(from).add(abIdx);
			outgoing.get(to).add(baIdx);
		}

		private static int nodeIndex(Coordinate c, List<Coordinate> nodes, Map<CoordKey, Integer> nodeLookup, List<List<Integer>> outgoing) {
			CoordKey key = new CoordKey(c.x, c.y);
			Integer idx = nodeLookup.get(key);
			if (idx != null) {
				Coordinate existing = nodes.get(idx.intValue());
				if (Double.isNaN(existing.getZ()) && !Double.isNaN(c.getZ())) {
					nodes.set(idx.intValue(), new Coordinate(existing.x, existing.y, c.getZ()));
				}
				return idx;
			}
			int newIndex = nodes.size();
			nodes.add(new Coordinate(c.x, c.y, c.getZ()));
			outgoing.add(new ArrayList<>());
			nodeLookup.put(key, newIndex);
			return newIndex;
		}

		private static Coordinate[] sanitizeRing(List<Coordinate> rawRing) {
			if (rawRing.size() < 3) {
				return null;
			}
			List<Coordinate> out = new ArrayList<>(rawRing.size() + 1);
			Coordinate prev = null;
			for (Coordinate c : rawRing) {
				if (prev == null || !prev.equals2D(c)) {
					out.add(c);
					prev = c;
				}
			}
			if (out.size() < 3) {
				return null;
			}
			Coordinate first = out.get(0);
			Coordinate last = out.get(out.size() - 1);
			if (!first.equals2D(last)) {
				out.add(new Coordinate(first.x, first.y, first.getZ()));
			}
			if (out.size() < 4) {
				return null;
			}
			return out.toArray(Coordinate[]::new);
		}

		private Geometry fallbackPolygonFaces(Polygon polygon) {
			Geometry coverage = polygon.getBoundary().union(asMultiLineString(polygon.getFactory()));
			var prepared = new IndexedPointInAreaLocator(polygon);

			Polygonizer polygonizer = new Polygonizer();
			polygonizer.setCheckRingsValid(false);
			polygonizer.add(coverage);
			List<Polygon> faces = new ArrayList<>();
			for (Object cell : polygonizer.getPolygons()) {
				if (cell instanceof Polygon face && prepared.locate(face.getInteriorPoint().getCoordinate()) == Location.INTERIOR) {
					faces.add(face);
				}
			}
			if (faces.isEmpty()) {
				return polygon.getFactory().createGeometryCollection();
			}
			return polygon.getFactory().buildGeometry(faces);
		}
	}
}
