package com.github.micycle1.grassfire4j.input;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.tinfour.common.IConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.HilbertSort;
import org.tinfour.utils.TriangleCollector;

import com.github.micycle1.grassfire4j.geom.Geom.Vec2;
import com.github.micycle1.grassfire4j.input.InputMesh.InputVertex;
import com.github.micycle1.grassfire4j.input.InputMeshBuilder.EdgeRef;

/**
 * Adapter that converts JTS {@link Polygon} input into {@link InputMesh}.
 * <p>
 * This is one concrete adapter for the solver boundary; additional adapters
 * for other input formats can be introduced by producing the same
 * {@link InputMesh} representation.
 */
public class PolygonAdapter implements Adapter<Polygon> {

	private record Edge(Vec2 a, Vec2 b) {
		Edge(Vec2 a, Vec2 b) {
			if (a.x() < b.x() || (a.x() == b.x() && a.y() < b.y())) {
				this.a = a;
				this.b = b;
			} else {
				this.a = b;
				this.b = a;
			}
		}
	}

	@Override
	public InputMesh toMesh(Polygon polygon) {
		return toMesh(polygon, null);
	}

	/**
	 * Converts a polygon into {@link InputMesh} while assigning per-boundary-edge
	 * weights.
	 * <p>
	 * Edge weights are consumed in boundary traversal order: exterior ring first,
	 * then interior rings in index order ({@code polygon.getInteriorRingN(i)}).
	 * Within each ring, each segment between consecutive coordinates contributes
	 * one weight entry (the closing duplicate coordinate is excluded).
	 *
	 * @param polygon input polygon
	 * @param edgeWeights optional per-edge weights; when non-null, size must equal
	 *        total boundary edge count across shell and holes
	 * @return triangulated solver input mesh
	 */
	public InputMesh toMesh(Polygon polygon, List<Double> edgeWeights) {
		if (polygon.isEmpty()) {
			return new InputMesh(List.of(), List.of());
		}

		List<Edge> boundaryEdges = getBoundaryEdgesInOrder(polygon);
		if (edgeWeights != null && edgeWeights.size() != boundaryEdges.size()) {
			throw new IllegalArgumentException(
					"edgeWeights size (" + edgeWeights.size() + ") must equal boundary edge count (" + boundaryEdges.size() + ")");
		}

		final var coords = polygon.getCoordinates();
		final double estPointSpacing = coords[0].distance(coords[1]);
		IncrementalTin tin = new IncrementalTin(estPointSpacing);

		List<Vertex> seedVertices = new ArrayList<>();
		Set<Vec2> uniqueVertices = new HashSet<>();
		addRingVertices(polygon.getExteriorRing(), seedVertices, uniqueVertices);
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			addRingVertices(polygon.getInteriorRingN(i), seedVertices, uniqueVertices);
		}

		if (seedVertices.size() < 3) {
			return new InputMesh(List.of(), List.of());
		}

		HilbertSort hilbertSort = new HilbertSort();
		hilbertSort.sort(seedVertices);
		tin.add(seedVertices, null);

		List<IConstraint> constraints = new ArrayList<>();
		PolygonConstraint exterior = createPolygonConstraint(polygon.getExteriorRing(), true);
		if (exterior != null) {
			constraints.add(exterior);
		}
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			PolygonConstraint hole = createPolygonConstraint(polygon.getInteriorRingN(i), false);
			if (hole != null) {
				constraints.add(hole);
			}
		}
		if (!constraints.isEmpty()) {
			tin.addConstraints(constraints, false);
		}

		Map<Vec2, Integer> vIndex = new HashMap<>();
		List<InputVertex> vertices = new ArrayList<>();
		List<int[]> triVertices = new ArrayList<>();

		TriangleCollector.visitTrianglesConstrained(tin, tri -> {
			int[] vIdx = new int[3];
			for (int i = 0; i < 3; i++) {
				Vec2 p = new Vec2(tri[i].getX(), tri[i].getY());
				vIdx[i] = vIndex.computeIfAbsent(p, pt -> {
					vertices.add(new InputVertex(pt.x(), pt.y(), true, null));
					return vertices.size() - 1;
				});
			}
			triVertices.add(vIdx);
		});

		Map<EdgeRef, Double> constrainedIndexEdgeWeights = new LinkedHashMap<>();
		for (int i = 0; i < boundaryEdges.size(); i++) {
			Edge edge = boundaryEdges.get(i);
			Integer a = vIndex.get(edge.a());
			Integer b = vIndex.get(edge.b());
			if (a != null && b != null) {
				EdgeRef edgeRef = new EdgeRef(a, b);
				double weight = edgeWeights == null ? 1.0 : edgeWeights.get(i);
				Double existing = constrainedIndexEdgeWeights.putIfAbsent(edgeRef, weight);
				if (existing != null && Double.compare(existing, weight) != 0) {
					throw new IllegalArgumentException("Conflicting weights mapped to the same boundary edge");
				}
			}
		}

		return InputMeshBuilder.build(vertices, triVertices, constrainedIndexEdgeWeights);
	}

	public static InputMesh fromPolygon(Polygon polygon) {
		return new PolygonAdapter().toMesh(polygon);
	}

	public static InputMesh fromPolygon(Polygon polygon, List<Double> edgeWeights) {
		return new PolygonAdapter().toMesh(polygon, edgeWeights);
	}

	private static List<Edge> getBoundaryEdgesInOrder(Polygon polygon) {
		List<Edge> boundaryEdges = new ArrayList<>();
		addRingEdges(polygon.getExteriorRing(), boundaryEdges);
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			addRingEdges(polygon.getInteriorRingN(i), boundaryEdges);
		}
		return boundaryEdges;
	}

	private static void addRingVertices(LineString ring, List<Vertex> out, Set<Vec2> seen) {
		Coordinate[] coords = ring.getCoordinates();
		for (int i = 0; i < coords.length - 1; i++) {
			Vec2 point = new Vec2(coords[i].x, coords[i].y);
			if (seen.add(point)) {
				out.add(new Vertex(point.x(), point.y(), Double.NaN));
			}
		}
	}

	private static PolygonConstraint createPolygonConstraint(LineString ring, boolean exterior) {
		Coordinate[] coords = ring.getCoordinates();
		if (coords.length < 4) {
			return null;
		}

		List<Vertex> points = new ArrayList<>(coords.length - 1);
		for (int i = 0; i < coords.length - 1; i++) {
			points.add(new Vertex(coords[i].x, coords[i].y, Double.NaN));
		}
		if (points.size() < 3) {
			return null;
		}

		boolean ccw = Orientation.isCCW(coords);
		if ((exterior && !ccw) || (!exterior && ccw)) {
			Collections.reverse(points);
		}

		return new PolygonConstraint(points);
	}

	private static void addRingEdges(LineString ring, List<Edge> edges) {
		Coordinate[] coords = ring.getCoordinates();
		for (int i = 0; i < coords.length - 1; i++) {
			edges.add(new Edge(new Vec2(coords[i].x, coords[i].y), new Vec2(coords[i+1].x, coords[i+1].y)));
		}
	}
}
