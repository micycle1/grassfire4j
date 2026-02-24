package com.github.micycle1.grassfire4j.input;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.micycle1.grassfire4j.input.InputMesh.Constraint;
import com.github.micycle1.grassfire4j.input.InputMesh.InputTriangle;
import com.github.micycle1.grassfire4j.input.InputMesh.InputVertex;

/**
 * Helper for building {@link InputMesh} from triangle vertex indices.
 * <p>
 * The builder computes triangle neighbor links and applies constrained edges,
 * so adapter implementations only need to provide vertices, triangles
 * (`int[3]`), and optional constrained edge index pairs.
 */
public final class InputMeshBuilder {

	private InputMeshBuilder() {
	}

	public static record EdgeRef(int a, int b) {
		public EdgeRef {
			if (a < 0 || b < 0) {
				throw new IllegalArgumentException("Edge indices must be non-negative");
			}
			if (a == b) {
				throw new IllegalArgumentException("Edge endpoints must be distinct");
			}
			if (a > b) {
				int t = a;
				a = b;
				b = t;
			}
		}
	}

	private static record SideRef(int tIdx, int side) {
	}

	public static InputMesh build(List<InputVertex> vertices, List<int[]> triangleVertices, Set<EdgeRef> constrainedEdges) {
		return build(vertices, triangleVertices, constrainedEdges, true);
	}

	public static InputMesh build(List<InputVertex> vertices, List<int[]> triangleVertices, Set<EdgeRef> constrainedEdges, boolean isInternal) {
		Map<EdgeRef, Double> edgeWeights = new HashMap<>();
		for (EdgeRef edgeRef : constrainedEdges) {
			edgeWeights.put(edgeRef, 1.0);
		}
		return build(vertices, triangleVertices, edgeWeights, isInternal);
	}

	public static InputMesh build(List<InputVertex> vertices, List<int[]> triangleVertices, Map<EdgeRef, Double> edgeWeights) {
		return build(vertices, triangleVertices, edgeWeights, true);
	}

	public static InputMesh build(List<InputVertex> vertices, List<int[]> triangleVertices, Map<EdgeRef, Double> edgeWeights, boolean isInternal) {
		Map<EdgeRef, Double> weights = edgeWeights == null ? Collections.emptyMap() : edgeWeights;
		Map<EdgeRef, List<SideRef>> edgeToSides = new HashMap<>();
		for (int tIdx = 0; tIdx < triangleVertices.size(); tIdx++) {
			int[] tv = triangleVertices.get(tIdx);
			validateTriangle(tv, vertices.size());
			for (int side = 0; side < 3; side++) {
				int a = tv[(side + 1) % 3];
				int b = tv[(side + 2) % 3];
				edgeToSides.computeIfAbsent(new EdgeRef(a, b), k -> new ArrayList<>()).add(new SideRef(tIdx, side));
			}
		}

		int[][] triN = new int[triangleVertices.size()][3];
		for (int[] row : triN) {
			Arrays.fill(row, -1);
		}
		for (List<SideRef> sides : edgeToSides.values()) {
			if (sides.size() == 2) {
				SideRef s0 = sides.get(0);
				SideRef s1 = sides.get(1);
				triN[s0.tIdx][s0.side] = s1.tIdx;
				triN[s1.tIdx][s1.side] = s0.tIdx;
			}
		}

		List<InputTriangle> triangles = new ArrayList<>(triangleVertices.size());
		for (int tIdx = 0; tIdx < triangleVertices.size(); tIdx++) {
			int[] tv = triangleVertices.get(tIdx);
			Constraint[] tc = new Constraint[3];
			for (int side = 0; side < 3; side++) {
				int a = tv[(side + 1) % 3];
				int b = tv[(side + 2) % 3];
				Double weight = weights.get(new EdgeRef(a, b));
				if (weight != null) {
					validateWeight(weight);
					tc[side] = new Constraint(weight);
				}
			}
			triangles.add(new InputTriangle(tv, triN[tIdx], tc, isInternal));
		}

		return new InputMesh(vertices, triangles);
	}

	private static void validateWeight(Double weight) {
		if (weight == null || !Double.isFinite(weight) || weight <= 0) {
			throw new IllegalArgumentException("Edge weight must be finite and > 0");
		}
	}

	private static void validateTriangle(int[] triangle, int vertexCount) {
		if (triangle == null || triangle.length != 3) {
			throw new IllegalArgumentException("Each triangle must be an int[3]");
		}
		int a = triangle[0];
		int b = triangle[1];
		int c = triangle[2];
		if (a < 0 || b < 0 || c < 0 || a >= vertexCount || b >= vertexCount || c >= vertexCount) {
			throw new IllegalArgumentException("Triangle contains out-of-range vertex index");
		}
		if (a == b || b == c || a == c) {
			throw new IllegalArgumentException("Triangle vertices must be distinct");
		}
	}
}
