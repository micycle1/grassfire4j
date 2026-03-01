package com.github.micycle1.grassfire4j.core;

import static com.github.micycle1.grassfire4j.geom.Geom.getBisector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.grassfire4j.geom.Geom.WaveFront;
import com.github.micycle1.grassfire4j.input.InputMesh;
import com.github.micycle1.grassfire4j.input.InputMesh.InputTriangle;
import com.github.micycle1.grassfire4j.input.InputMesh.InputVertex;
import com.github.micycle1.grassfire4j.model.Model.InfiniteVertex;
import com.github.micycle1.grassfire4j.model.Model.KineticTriangle;
import com.github.micycle1.grassfire4j.model.Model.KineticVertex;
import com.github.micycle1.grassfire4j.model.Model.Skeleton;
import com.github.micycle1.grassfire4j.model.Model.SkeletonNode;

/**
 * Core initialisation routines that convert an {@link InputMesh} into the
 * kinetic model used by the event-driven straight-skeleton solver.
 * <p>
 * Since this layer depends only on {@link InputMesh}, it is decoupled from any
 * particular frontend adapter.
 */
public final class Core {

	private Core() {
	}

	public static int ccw(int i) { return (i + 1) % 3; }
	public static int cw(int i) { return (i + 2) % 3; }
	public static int apex(int i) { return i % 3; }
	public static int orig(int i) { return (i + 1) % 3; }
	public static int dest(int i) { return (i + 2) % 3; }

	public record Corner(int tIdx, int side) {}

	private record PointKey(double x, double y) implements Comparable<PointKey> {
		@Override
		public int compareTo(PointKey o) {
			int cx = Double.compare(this.x, o.x);
			return cx != 0 ? cx : Double.compare(this.y, o.y);
		}
	}

	private record TriangleKey(PointKey a, PointKey b, PointKey c) implements Comparable<TriangleKey> {
		static TriangleKey from(InputMesh mesh, InputTriangle tri) {
			PointKey[] pts = new PointKey[3];
			for (int i = 0; i < 3; i++) {
				InputVertex v = mesh.vertices.get(tri.v[i]);
				pts[i] = new PointKey(v.x, v.y);
			}
			Arrays.sort(pts);
			return new TriangleKey(pts[0], pts[1], pts[2]);
		}

		@Override
		public int compareTo(TriangleKey o) {
			int ca = this.a.compareTo(o.a);
			if (ca != 0) {
				return ca;
			}
			int cb = this.b.compareTo(o.b);
			if (cb != 0) {
				return cb;
			}
			return this.c.compareTo(o.c);
		}
	}

	private static int cornerId(int tIdx, int side) { return tIdx * 3 + side; }
	private static int cornerTri(int cornerId) { return cornerId / 3; }
	private static int cornerSide(int cornerId) { return cornerId % 3; }

	private static int stepAroundVertex(InputMesh mesh, int vIdx, int cornerId, boolean isCcw) {
		int tIdx = cornerTri(cornerId);
		int side = cornerSide(cornerId);
		InputTriangle tri = mesh.triangles.get(tIdx);
		int edgeSide = isCcw ? ccw(side) : cw(side);
		if (tri.c[edgeSide] != null) {
			return -1;
		}
		int nIdx = tri.n[edgeSide];
		if (nIdx < 0) {
			return -1;
		}
		int[] vArr = mesh.triangles.get(nIdx).v;
		int nSide = -1;
		for (int i = 0; i < 3; i++) {
			if (vArr[i] == vIdx) {
				nSide = i;
				break;
			}
		}
		return nSide < 0 ? -1 : cornerId(nIdx, nSide);
	}

	public static List<List<List<Corner>>> buildVertexStars(InputMesh mesh) {
		int vertexCount = mesh.vertices.size();
		int triangleCount = mesh.triangles.size();
		int totalCorners = triangleCount * 3;

		List<List<Integer>> stars = new ArrayList<>(vertexCount);
		List<List<List<Corner>>> grouped = new ArrayList<>(vertexCount);
		for (int i = 0; i < vertexCount; i++) {
			stars.add(new ArrayList<>());
			grouped.add(new ArrayList<>());
		}
		for (int tIdx = 0; tIdx < triangleCount; tIdx++) {
			int[] vArr = mesh.triangles.get(tIdx).v;
			for (int side = 0; side < 3; side++) {
				stars.get(vArr[side]).add(cornerId(tIdx, side));
			}
		}

		int[] unvisitedMark = new int[totalCorners];
		int[] seenMark = new int[totalCorners];
		int[] walkedMark = new int[totalCorners];
		int epoch = 1;

		for (int vIdx = 0; vIdx < stars.size(); vIdx++) {
			List<Integer> incidentCorners = stars.get(vIdx);
			int unvisitedEpoch = epoch++;
			int remaining = incidentCorners.size();
			for (int c : incidentCorners) {
				unvisitedMark[c] = unvisitedEpoch;
			}

			while (remaining > 0) {
				int start = -1;
				for (int c : incidentCorners) {
					if (unvisitedMark[c] == unvisitedEpoch) {
						start = c;
						break;
					}
				}
				if (start < 0) {
					break;
				}

				int seenEpoch = epoch++;
				int leftmost = start;
				seenMark[start] = seenEpoch;
				while (true) {
					int prev = stepAroundVertex(mesh, vIdx, leftmost, true);
					if (prev < 0 || prev == leftmost || seenMark[prev] == seenEpoch) {
						break;
					}
					seenMark[prev] = seenEpoch;
					leftmost = prev;
				}

				List<Corner> group = new ArrayList<>();
				int walkedEpoch = epoch++;
				int cur = leftmost;
				while (cur >= 0 && walkedMark[cur] != walkedEpoch) {
					walkedMark[cur] = walkedEpoch;
					group.add(new Corner(cornerTri(cur), cornerSide(cur)));
					if (unvisitedMark[cur] == unvisitedEpoch) {
						unvisitedMark[cur] = 0;
						remaining--;
					}
					int next = stepAroundVertex(mesh, vIdx, cur, false);
					if (leftmost == next) {
						break;
					}
					cur = next;
				}
				if (!group.isEmpty()) {
					Collections.reverse(group);
					grouped.get(vIdx).add(group);
				}
			}
		}
		return grouped;
	}

	public static Skeleton initSkeleton(InputMesh mesh) {
		Skeleton skel = new Skeleton();
		Map<Integer, SkeletonNode> nodes = new HashMap<>();
		Map<Integer, InfiniteVertex> infNodes = new HashMap<>();
		Map<TriangleKey, Integer> canonicalUids = canonicalTriangleUids(mesh);
		double sumX = 0, sumY = 0;
		int count = 0;

		for (int i = 0; i < mesh.vertices.size(); i++) {
			InputVertex v = mesh.vertices.get(i);
			if (v.isFinite) {
				nodes.put(i, new SkeletonNode(new Vector2D(v.x, v.y), -1, v.info));
				sumX += v.x;
				sumY += v.y;
				count++;
			} else {
				infNodes.put(i, new InfiniteVertex(new Vector2D(v.x, v.y)));
			}
		}
		InfiniteVertex centroid = new InfiniteVertex(new Vector2D(count > 0 ? sumX / count : 0, count > 0 ? sumY / count : 0));

		List<KineticTriangle> ktriangles = new ArrayList<>();
		for (int i = 0; i < mesh.triangles.size(); i++) {
			KineticTriangle k = new KineticTriangle();
			k.info = i + 1;
			k.uid = canonicalUids.get(TriangleKey.from(mesh, mesh.triangles.get(i))).intValue();
			k.internal = mesh.triangles.get(i).isInternal;
			ktriangles.add(k);
		}

		for (int tIdx = 0; tIdx < mesh.triangles.size(); tIdx++) {
			InputTriangle t = mesh.triangles.get(tIdx);
			KineticTriangle k = ktriangles.get(tIdx);
			for (int i = 0; i < 3; i++) {
				if (!mesh.vertices.get(t.v[i]).isFinite) {
					k.vertices[i] = infNodes.get(t.v[i]);
				}
			}
			for (int i = 0; i < 3; i++) {
				if (t.c[i] != null) {
					InputVertex start = mesh.vertices.get(t.v[ccw(i)]);
					InputVertex end = mesh.vertices.get(t.v[cw(i)]);
					Object wavefrontData = t.c[i].edgeId() >= 0 ? t.c[i].edgeId() : null;
					k.wavefrontSupportLines[i] = new WaveFront(new Vector2D(start.x, start.y), new Vector2D(end.x, end.y), null, t.c[i].weight(),
							wavefrontData);
				}
			}
			for (int i = 0; i < 3; i++) {
				if (t.c[i] == null && t.n[i] != -1) {
					k.neighbours[i] = ktriangles.get(t.n[i]);
				}
			}
		}

		var stars = buildVertexStars(mesh);
		int ct = 0;
		record Link(Corner left, KineticVertex kv, Corner right) {
		}
		List<Link> links = new ArrayList<>();

		for (int vIdx = 0; vIdx < stars.size(); vIdx++) {
			if (!mesh.vertices.get(vIdx).isFinite) {
				continue;
			}
			for (List<Corner> group : stars.get(vIdx)) {
				Corner first = group.get(0), last = group.get(group.size() - 1);
				InputTriangle tLast = mesh.triangles.get(last.tIdx()), tFirst = mesh.triangles.get(first.tIdx());
				InputVertex tail = mesh.vertices.get(tLast.v[cw(last.side())]);
				InputVertex mid = mesh.vertices.get(tLast.v[last.side()]);
				InputVertex head = mesh.vertices.get(tFirst.v[ccw(first.side())]);

				int turn = Orientation.index(new Coordinate(tail.x, tail.y), new Coordinate(mid.x, mid.y), new Coordinate(head.x, head.y));
				KineticVertex.Turn turnType = turn == -1 ? KineticVertex.Turn.RIGHT_REFLEX
						: (turn == 1 ? KineticVertex.Turn.LEFT_CONVEX : KineticVertex.Turn.STRAIGHT);

				WaveFront right = ktriangles.get(first.tIdx()).wavefrontSupportLines[cw(first.side())];
				WaveFront left = ktriangles.get(last.tIdx()).wavefrontSupportLines[ccw(last.side())];

				KineticVertex kv = new KineticVertex();
				kv.info = ++ct;
				kv.origin = new Vector2D(mesh.vertices.get(vIdx).x, mesh.vertices.get(vIdx).y);
				kv.velocity = getBisector(left, right);
				kv.startNode = nodes.get(vIdx);
				kv.startsAt = 0.0;
				kv.turn = turnType;
				kv.ul = left.line;
				kv.ur = right.line;
				kv.wfl = left;
				kv.wfr = right;

				for (Corner c : group) {
					KineticTriangle kt = ktriangles.get(c.tIdx());
					kt.vertices[c.side()] = kv;
					kv.internal = kt.internal;
				}
				skel.vertices.add(kv);
				links.add(new Link(new Corner(last.tIdx(), cw(last.side())), kv, new Corner(first.tIdx(), ccw(first.side()))));
			}
		}

		for (Link l : links) {
			KineticVertex lkv = (KineticVertex) ktriangles.get(l.left().tIdx()).vertices[l.left().side()];
			KineticVertex rkv = (KineticVertex) ktriangles.get(l.right().tIdx()).vertices[l.right().side()];
			l.kv().setLeft(lkv, 0.0);
			l.kv().setRight(rkv, 0.0);
		}

		for (KineticTriangle t : ktriangles) {
			for (int i = 0; i < 3; i++) {
				if (t.vertices[i] instanceof InfiniteVertex) {
					t.vertices[i] = centroid;
				}
			}
		}

		ktriangles.sort(Comparator.comparingDouble((KineticTriangle t) -> t.vertices[0].positionAt(0).getY())
				.thenComparingDouble(t -> t.vertices[0].positionAt(0).getX()));
		skel.skNodes.addAll(nodes.values());
		skel.triangles = ktriangles;
		skel.boundaryEdges = buildBoundaryEdgesFromMesh(mesh);
		return skel;
	}

	private static Map<TriangleKey, Integer> canonicalTriangleUids(InputMesh mesh) {
		Map<TriangleKey, Integer> out = new HashMap<>();
		TreeMap<TriangleKey, Integer> ordered = new TreeMap<>();
		for (InputTriangle tri : mesh.triangles) {
			ordered.putIfAbsent(TriangleKey.from(mesh, tri), Integer.valueOf(0));
		}
		int uid = 1;
		for (TriangleKey key : ordered.keySet()) {
			out.put(key, Integer.valueOf(uid++));
		}
		return out;
	}

	/**
	 * Extracts directed boundary edges from the triangulated input mesh.
	 * <p>
	 * Edges are oriented consistently such that the interior of the polygon
	 * (represented by internal triangles) always lies to the <i>left</i> of the
	 * directed edge (pointing from {@code from} to {@code to}).
	 * <p>
	 * The output list is returned in ascending order of their assigned
	 * {@code edgeId}. This ensures a stable boundary representation that
	 * corresponds predictably to the original input segment ordering.
	 *
	 * @param mesh the populated input mesh containing vertices and triangles
	 * @return a deterministically sorted list of oriented boundary edges
	 */
	private static List<Skeleton.BoundaryEdge> buildBoundaryEdgesFromMesh(InputMesh mesh) {
		Map<Integer, Skeleton.BoundaryEdge> byId = new HashMap<>();

		for (InputTriangle t : mesh.triangles) {
			if (!t.isInternal) { // important: orient using an interior triangle
				continue;
			}

			for (int i = 0; i < 3; i++) {
				if (t.c[i] == null) {
					continue;
				}
				int edgeId = t.c[i].edgeId();
				if (edgeId < 0) {
					continue;
				}

				int aIdx = t.v[ccw(i)];
				int bIdx = t.v[cw(i)];
				int cIdx = t.v[i]; // third vertex of the triangle

				InputVertex av = mesh.vertices.get(aIdx);
				InputVertex bv = mesh.vertices.get(bIdx);
				InputVertex cv = mesh.vertices.get(cIdx);

				Coordinate a = new Coordinate(av.x, av.y);
				Coordinate b = new Coordinate(bv.x, bv.y);
				Coordinate c = new Coordinate(cv.x, cv.y);

				// Orient boundary edge so the interior triangle lies on the LEFT of (from ->
				// to).
				int turn = Orientation.index(a, b, c);
				if (turn == 0) {
					continue;
				}

				Coordinate from = (turn > 0) ? a : b;
				Coordinate to = (turn > 0) ? b : a;

				byId.putIfAbsent(edgeId, new Skeleton.BoundaryEdge(from, to, edgeId));
			}
		}

		// return in edgeId order (stable)
		if (byId.isEmpty()) {
			return List.of();
		}
		int max = byId.keySet().stream().mapToInt(Integer::intValue).max().orElse(-1);

		List<Skeleton.BoundaryEdge> out = new ArrayList<>(byId.size());
		for (int id = 0; id <= max; id++) {
			Skeleton.BoundaryEdge e = byId.get(id);
			if (e != null) {
				out.add(e);
			}
		}
		return out;
	}
}
