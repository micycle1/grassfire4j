package grassfire2;

import grassfire2.core.*;
import grassfire2.events.Events;
import grassfire2.model.Model.Skeleton;
import grassfire2.triangulation.JtsAdapter;
import grassfire2.triangulation.JtsAdapter.InputMesh;

import org.locationtech.jts.geom.Polygon;

public class GrassfireApi {

	public static Skeleton computeSkeleton(Polygon polygon) {
		// 1. Triangulate Polygon into Internal Mesh
		var mesh = JtsAdapter.fromJtsPolygon(polygon);
		return computeSkeleton(mesh);
	}

	public static Skeleton computeSkeletonTinfour(Polygon polygon) {
		// 1. Triangulate Polygon into Internal Mesh
		var mesh = JtsAdapter.fromJtsPolygonTinfour(polygon);
		return computeSkeleton(mesh);
	}

	public static Skeleton computeSkeleton(InputMesh mesh) {
		// 2. Init Kinetic Skeleton Data Structure
		Skeleton skel = Core.initSkeleton(mesh);

		// 3. Initialize Event Queue
		Events.EventQueue queue = new Events.EventQueue();
		for (var tri : skel.triangles) {
			var e = CollapseEventComputer.compute(tri, 0.0, true);
			if (e != null)
				queue.add(e);
		}

		// 4. Run Physics Loop
		Events.eventLoop(queue, skel);

		return skel;
	}

}
