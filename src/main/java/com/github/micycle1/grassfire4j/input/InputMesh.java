package com.github.micycle1.grassfire4j.input;

import java.util.List;

/**
 * Intermediate mesh representation consumed by the kinetic skeleton solver.
 * <p>
 * Input adapters should translate source geometry into this model
 * ({@link InputVertex}, {@link InputTriangle}, and optional edge
 * {@link Constraint} data). The solver core operates on this type only, which
 * makes it straightforward to add new adapters without changing solver logic.
 */
public class InputMesh {

	/**
	 * Per-edge constraint metadata.
	 *
	 * @param weight wavefront speed weight for the constrained edge
	 * @param edgeId optional stable boundary edge identifier ({@code -1} if absent)
	 */
	public record Constraint(double weight, int edgeId) {
		public Constraint(double weight) {
			this(weight, -1);
		}
	}

	public static class InputVertex {
		public final double x, y;
		public final boolean isFinite;
		public final Integer info;

		public InputVertex(double x, double y, boolean isFinite, Integer info) {
			this.x = x;
			this.y = y;
			this.isFinite = isFinite;
			this.info = info;
		}
	}

	public static class InputTriangle {
		public final int[] v;
		public final int[] n;
		public final Constraint[] c;
		public final boolean isInternal;

		public InputTriangle(int[] v, int[] n, Constraint[] c, boolean isInternal) {
			this.v = v;
			this.n = n;
			this.c = c;
			this.isInternal = isInternal;
		}
	}

	public final List<InputVertex> vertices;
	public final List<InputTriangle> triangles;

	public InputMesh(List<InputVertex> vertices, List<InputTriangle> triangles) {
		this.vertices = vertices;
		this.triangles = triangles;
	}
}
