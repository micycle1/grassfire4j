package com.github.micycle1.grassfire4j.input;

/**
 * Converts an external input type into the solver's intermediate
 * {@link InputMesh} representation.
 * <p>
 * {@link #toMesh(Object)} must return triangles that are fully wired for the
 * solver:
 * </p>
 * <p>
 * <b>Triangle connectivity contract:</b>
 * </p>
 * <ul>
 * <li>{@code InputTriangle.v} stores 3 vertex indices.</li>
 * <li>{@code InputTriangle.n} stores 3 neighbor triangle indices ({@code -1}
 * for boundary).</li>
 * <li>{@code InputTriangle.c} stores 3 edge constraints ({@code null} for
 * unconstrained).</li>
 * <li>Triangle side {@code i} is opposite vertex {@code v[i]}; its edge
 * endpoints are {@code v[(i+1)%3]} and {@code v[(i+2)%3]}.</li>
 * </ul>
 * <p>
 * Implementations may construct these arrays directly, or use
 * {@link InputMeshBuilder} to derive neighbor links and constrained edges from
 * triangle vertex indices.
 * </p>
 *
 * @param <T> source input type
 */
public interface Adapter<T> {

	InputMesh toMesh(T input);
}
