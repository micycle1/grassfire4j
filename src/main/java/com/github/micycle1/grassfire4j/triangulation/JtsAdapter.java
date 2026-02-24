package com.github.micycle1.grassfire4j.triangulation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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

public class JtsAdapter {

	public record Constraint(double weight) {}

	public static class InputVertex {
		public final double x, y;
		public final boolean isFinite;
		public final Integer info;
		public InputVertex(double x, double y, boolean isFinite, Integer info) {
			this.x = x; this.y = y; this.isFinite = isFinite; this.info = info;
		}
	}

	public static class InputTriangle {
		public final int[] v;
		public final int[] n;
		public final Constraint[] c;
		public final boolean isInternal;
		public InputTriangle(int[] v, int[] n, Constraint[] c, boolean isInternal) {
			this.v = v; this.n = n; this.c = c; this.isInternal = isInternal;
		}
	}

	public static class InputMesh {
		public final List<InputVertex> vertices;
		public final List<InputTriangle> triangles;
		public InputMesh(List<InputVertex> vertices, List<InputTriangle> triangles) {
			this.vertices = vertices; this.triangles = triangles;
		}
	}

	private record Edge(Vec2 a, Vec2 b) {
		Edge(Vec2 a, Vec2 b) {
			if (a.x() < b.x() || (a.x() == b.x() && a.y() < b.y())) { this.a = a; this.b = b; }
			else { this.a = b; this.b = a; }
		}
	}

	public static InputMesh fromPolygon(Polygon polygon) {
		if (polygon.isEmpty()) {
			return new InputMesh(List.of(), List.of());
		}

		Set<Edge> constrainedEdges = getConstrainedEdges(polygon);

		IncrementalTin tin = new IncrementalTin(1.0);
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

		return buildInputMesh(vertices, triVertices, constrainedEdges);
	}

	private static InputMesh buildInputMesh(List<InputVertex> vertices, List<int[]> triVertices, Set<Edge> constrainedEdges) {
		record SideRef(int tIdx, int side) {}
		Map<Edge, List<SideRef>> edgeToSides = new HashMap<>();
		for (int tIdx = 0; tIdx < triVertices.size(); tIdx++) {
			int[] tv = triVertices.get(tIdx);
			for (int side = 0; side < 3; side++) {
				InputVertex a = vertices.get(tv[(side + 1) % 3]);
				InputVertex b = vertices.get(tv[(side + 2) % 3]); // side-1 mapped to side+2
				edgeToSides.computeIfAbsent(new Edge(new Vec2(a.x, a.y), new Vec2(b.x, b.y)), k -> new ArrayList<>()).add(new SideRef(tIdx, side));
			}
		}

		List<InputTriangle> triangles = new ArrayList<>();
		int[][] triN = new int[triVertices.size()][3];
		for (var row : triN) {
			Arrays.fill(row, -1);
		}

		for (var sides : edgeToSides.values()) {
			if (sides.size() == 2) {
				SideRef s0 = sides.get(0), s1 = sides.get(1);
				triN[s0.tIdx][s0.side] = s1.tIdx;
				triN[s1.tIdx][s1.side] = s0.tIdx;
			}
		}

		for (int tIdx = 0; tIdx < triVertices.size(); tIdx++) {
			int[] tv = triVertices.get(tIdx);
			Constraint[] tc = new Constraint[3];
			for (int side = 0; side < 3; side++) {
				InputVertex a = vertices.get(tv[(side + 1) % 3]);
				InputVertex b = vertices.get(tv[(side + 2) % 3]);
				if (constrainedEdges.contains(new Edge(new Vec2(a.x, a.y), new Vec2(b.x, b.y)))) {
					tc[side] = new Constraint(1.0);
				}
			}
			triangles.add(new InputTriangle(tv, triN[tIdx], tc, true));
		}

		return new InputMesh(vertices, triangles);
	}

	private static Set<Edge> getConstrainedEdges(Polygon polygon) {
		Set<Edge> constrainedEdges = new HashSet<>();
		addRingEdges(polygon.getExteriorRing(), constrainedEdges);
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			addRingEdges(polygon.getInteriorRingN(i), constrainedEdges);
		}
		return constrainedEdges;
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

	private static void addRingEdges(LineString ring, Set<Edge> edges) {
		Coordinate[] coords = ring.getCoordinates();
		for (int i = 0; i < coords.length - 1; i++) {
			edges.add(new Edge(new Vec2(coords[i].x, coords[i].y), new Vec2(coords[i+1].x, coords[i+1].y)));
		}
	}
}
