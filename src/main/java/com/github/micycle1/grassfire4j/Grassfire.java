package com.github.micycle1.grassfire4j;

import org.locationtech.jts.geom.Polygon;

import com.github.micycle1.grassfire4j.core.CollapseEventComputer;
import com.github.micycle1.grassfire4j.core.Core;
import com.github.micycle1.grassfire4j.events.Events;
import com.github.micycle1.grassfire4j.model.Model.Skeleton;
import com.github.micycle1.grassfire4j.triangulation.JtsAdapter;
import com.github.micycle1.grassfire4j.triangulation.JtsAdapter.InputMesh;

public class Grassfire {

	public static Skeleton computeSkeleton(Polygon polygon) {
		var mesh = JtsAdapter.fromPolygon(polygon);
		return computeSkeleton(mesh);
	}

	public static Skeleton computeSkeletonTinfour(Polygon polygon) {
		var mesh = JtsAdapter.fromPolygon(polygon);
		return computeSkeleton(mesh);
	}

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
