package com.github.micycle1.grassfire4j;

import java.util.Objects;

import org.locationtech.jts.geom.Polygon;

import com.github.micycle1.grassfire4j.core.CollapseEventComputer;
import com.github.micycle1.grassfire4j.core.Core;
import com.github.micycle1.grassfire4j.events.Events;
import com.github.micycle1.grassfire4j.input.Adapter;
import com.github.micycle1.grassfire4j.input.InputMesh;
import com.github.micycle1.grassfire4j.input.PolygonAdapter;
import com.github.micycle1.grassfire4j.model.Model.Skeleton;

/**
 * Public API for computing a kinetic straight skeleton (wavefront-collapse /
 * grassfire).
 * <p>
 * The API is intentionally split into adapters and solver core: adapters map
 * external geometry formats into {@link InputMesh}, and the skeleton solver
 * then consumes that intermediate representation.
 * <p>
 * This design allows additional adapters to be added without changing the
 * kinetic skeleton pipeline.
 */
public class Grassfire {

	private Grassfire() {
	}

	/**
	 * Computes the kinetic straight skeleton of a JTS {@link Polygon}.
	 * <p>
	 * Supports polygons with holes.
	 *
	 * @param polygon polygon input (outer shell and optional interior rings)
	 * @return skeleton model containing nodes, kinetic vertices, and segments
	 */
	public static Skeleton computeSkeleton(Polygon polygon) {
		return computeSkeleton(polygon, new PolygonAdapter());
	}

	/**
	 * Computes the kinetic straight skeleton from user-supplied input via an
	 * adapter.
	 *
	 * @param <T>     source input type
	 * @param input   source input object
	 * @param adapter input adapter that converts {@code input} into
	 *                {@link InputMesh}
	 * @return skeleton model containing nodes, kinetic vertices, and segments
	 */
	public static <T> Skeleton computeSkeleton(T input, Adapter<T> adapter) {
		Objects.requireNonNull(adapter, "adapter");
		return computeSkeleton(adapter.toMesh(input));
	}

	/**
	 * Computes the kinetic straight skeleton from a pre-built input mesh.
	 * <p>
	 * Use this overload when you need control over triangulation and per-edge
	 * weights before running the event-based solver.
	 *
	 * @param mesh triangulation-ready input mesh for the kinetic solver
	 * @return skeleton model containing nodes, kinetic vertices, and segments
	 */
	public static Skeleton computeSkeleton(InputMesh mesh) {
		Skeleton skel = Core.initSkeleton(mesh);

		Events.EventQueue queue = new Events.EventQueue();
		for (var tri : skel.triangles) {
			var e = CollapseEventComputer.compute(tri, 0.0, true);
			if (e != null) {
				queue.add(e);
			}
		}

		Events.eventLoop(queue, skel);

		return skel;
	}

}
