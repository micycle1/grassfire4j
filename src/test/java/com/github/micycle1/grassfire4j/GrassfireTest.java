package com.github.micycle1.grassfire4j;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;
import static org.junit.jupiter.api.DynamicTest.dynamicTest;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;
import java.util.stream.Collectors;

import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestFactory;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.coverage.CoverageValidator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.index.hprtree.HPRtree;

import com.github.micycle1.grassfire4j.input.InputMesh;
import com.github.micycle1.grassfire4j.input.PolygonAdapter;
import com.github.micycle1.grassfire4j.model.Model.Skeleton.Segment;

class GrassfireTest {

	private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();
	private static final Path CSV_DIR = Path.of("src", "test", "resources", "csv");
	private static final Map<String, Integer> EXPECTED_SEGMENTS = expectedSegments();

	@Test
	void internalSegmentsCount() {
		List<List<Coordinate>> rings = List.of(List.of(c(0.0, 0.0), c(20.0, 0.0), c(20.0, 10.0), c(10.0, 10.0), c(10.0, 20.0), c(0.0, 20.0), c(0.0, 0.0)));

		Polygon polygon = toPolygon(rings);
		var skeleton = Grassfire.computeSkeleton(polygon);
		assertEquals(8, skeleton.segments().size(), "Unexpected internal segment count");
	}

	@TestFactory
	Stream<DynamicTest> skeletonIntegrity() throws IOException {
		return csvFiles().map(csvFile -> dynamicTest(csvFile.getFileName().toString(), () -> runSkeletonIntegrity(csvFile)));
	}

	@Test
	void weightedPolygonDiffersFromUnweighted() {
		List<List<Coordinate>> rings = List.of(List.of(c(0.0, 0.0), c(10.0, 0.0), c(10.0, 6.0), c(0.0, 6.0), c(0.0, 0.0)));
		Polygon polygon = toPolygon(rings);

		var unweighted = Grassfire.computeSkeleton(polygon);
		var weighted = Grassfire.computeSkeleton(polygon, p -> new PolygonAdapter().toMesh(p, List.of(0.5, 3.0, 1.0, 1.0)));

		Set<String> unweightedSegments = canonicalSegments(unweighted.segments());
		Set<String> weightedSegments = canonicalSegments(weighted.segments());

		assertFalse(unweightedSegments.equals(weightedSegments), "Weighted and unweighted skeletons should differ");
	}

	@Test
	void boundaryVertexInfoAutoAssignedByIndex() {
		List<List<Coordinate>> rings = List.of(List.of(c(0.0, 0.0), c(10.0, 0.0), c(10.0, 10.0), c(0.0, 10.0), c(0.0, 0.0)));
		Polygon polygon = toPolygon(rings);

		var mesh = new PolygonAdapter().toMesh(polygon);

		assertEquals(Integer.valueOf(0), findVertexInfo(mesh, 0.0, 0.0));
		assertEquals(Integer.valueOf(1), findVertexInfo(mesh, 10.0, 0.0));
		assertEquals(Integer.valueOf(2), findVertexInfo(mesh, 10.0, 10.0));
		assertEquals(Integer.valueOf(3), findVertexInfo(mesh, 0.0, 10.0));
	}

	@Test
	void polygonFacesPreserveInputHole() {
		List<List<Coordinate>> rings = List.of(List.of(c(0.0, 0.0), c(12.0, 0.0), c(12.0, 12.0), c(0.0, 12.0), c(0.0, 0.0)),
				List.of(c(4.0, 4.0), c(8.0, 4.0), c(8.0, 8.0), c(4.0, 8.0), c(4.0, 4.0)));

		Polygon polygon = toPolygon(rings);
		var skeleton = Grassfire.computeSkeleton(polygon);
		var faces = skeleton.asPolygonFaces();

		assertEquals(polygon.getArea(), faces.getArea(), 1e-9, "Faces area should match polygon area (holes preserved)");
		assertFalse(faces.covers(GEOMETRY_FACTORY.createPoint(c(6.0, 6.0))), "Hole interior should remain empty");
	}

	@Test
	void polygonDifferentOrientations() throws IOException {
		Path csvPath = CSV_DIR.resolve("matisse-alga.csv");
		List<List<Coordinate>> rings = readCsvPolygon(csvPath);

		List<Coordinate> shell = rings.get(0);
		List<List<Coordinate>> holes = rings.subList(1, rings.size());

		List<Coordinate> shellCW = new ArrayList<>(shell);
		if (Orientation.isCCW(shellCW.toArray(new Coordinate[0]))) {
			Collections.reverse(shellCW);
		}
		List<Coordinate> shellCCW = new ArrayList<>(shell);
		if (!Orientation.isCCW(shellCCW.toArray(new Coordinate[0]))) {
			Collections.reverse(shellCCW);
		}

		List<List<Coordinate>> holesCW = new ArrayList<>();
		for (List<Coordinate> hole : holes) {
			List<Coordinate> holeCW = new ArrayList<>(hole);
			if (Orientation.isCCW(holeCW.toArray(new Coordinate[0]))) {
				Collections.reverse(holeCW);
			}
			holesCW.add(holeCW);
		}

		List<List<Coordinate>> holesCCW = new ArrayList<>();
		for (List<Coordinate> hole : holes) {
			List<Coordinate> holeCCW = new ArrayList<>(hole);
			if (!Orientation.isCCW(holeCCW.toArray(new Coordinate[0]))) {
				Collections.reverse(holeCCW);
			}
			holesCCW.add(holeCCW);
		}

		// CCW shell, CW holes
		List<List<Coordinate>> rings1 = new ArrayList<>();
		rings1.add(shellCCW);
		rings1.addAll(holesCW);
		Polygon polygon1 = toPolygon(rings1);

		// CW shell, CCW holes
		List<List<Coordinate>> rings2 = new ArrayList<>();
		rings2.add(shellCW);
		rings2.addAll(holesCCW);
		Polygon polygon2 = toPolygon(rings2);

		// CCW shell, CCW holes
		List<List<Coordinate>> rings3 = new ArrayList<>();
		rings3.add(shellCCW);
		rings3.addAll(holesCCW);
		Polygon polygon3 = toPolygon(rings3);

		// CW shell, CW holes
		List<List<Coordinate>> rings4 = new ArrayList<>();
		rings4.add(shellCW);
		rings4.addAll(holesCW);
		Polygon polygon4 = toPolygon(rings4);

		var skeleton1 = Grassfire.computeSkeleton(polygon1);
		var skeleton2 = Grassfire.computeSkeleton(polygon2);
		var skeleton3 = Grassfire.computeSkeleton(polygon3);
		var skeleton4 = Grassfire.computeSkeleton(polygon4);

		assertEquals(skeleton1.segments().size(), skeleton2.segments().size(), "Skeleton segment counts should match regardless of ring orientation");
		assertEquals(skeleton1.segments().size(), skeleton3.segments().size(), "Skeleton segment counts should match regardless of ring orientation");
		assertEquals(skeleton1.segments().size(), skeleton4.segments().size(), "Skeleton segment counts should match regardless of ring orientation");

		Set<String> segments1 = canonicalSegments(skeleton1.segments());
		Set<String> segments2 = canonicalSegments(skeleton2.segments());
		Set<String> segments3 = canonicalSegments(skeleton3.segments());
		Set<String> segments4 = canonicalSegments(skeleton4.segments());

		assertEquals(segments1, segments2, "Skeletons should be identical regardless of ring orientation");
		assertEquals(segments1, segments3, "Skeletons should be identical regardless of ring orientation");
		assertEquals(segments1, segments4, "Skeletons should be identical regardless of ring orientation");

		var sp1 = skeleton1.asPolygonFaces();
		var sp2 = skeleton2.asPolygonFaces();
		var sp3 = skeleton3.asPolygonFaces();
		var sp4 = skeleton4.asPolygonFaces();
		assertTrue(sp1.getNumPoints() == sp2.getNumPoints());
		assertTrue(sp3.getNumPoints() == sp4.getNumPoints());
		assertTrue(sp1.getNumPoints() == sp3.getNumPoints());
		assertTrue(sp1.getNumGeometries() == sp2.getNumGeometries());
		assertTrue(sp3.getNumGeometries() == sp4.getNumGeometries());
		assertTrue(sp1.getNumGeometries() == sp3.getNumGeometries());
		assertEquals(sp1.getArea(), sp2.getArea(), 1e-9);
		assertEquals(sp3.getArea(), sp4.getArea(), 1e-9);
		assertEquals(sp1.getArea(), sp3.getArea(), 1e-9);
	}

	private void runSkeletonIntegrity(Path csvFile) throws IOException {
		List<List<Coordinate>> rings = readCsvPolygon(csvFile);
		if (rings.isEmpty()) {
			return;
		}

		Set<PointKey> inputVertices = new HashSet<>();
		for (List<Coordinate> ring : rings) {
			for (Coordinate p : ring) {
				inputVertices.add(new PointKey(p.x, p.y));
			}
		}

		Polygon polygon = toPolygon(rings);
		var skeleton = Grassfire.computeSkeleton(polygon);
		List<Segment> segments = skeleton.segments();

		Integer expected = EXPECTED_SEGMENTS.get(csvFile.getFileName().toString());
		if (expected != null) {
			assertEquals(expected.intValue(), segments.size(), "Unexpected segment count for " + csvFile.getFileName());
		}

		Set<PointKey> skeletonEndpoints = new HashSet<>();
		List<SegmentLine> skeletonLines = new ArrayList<>(segments.size());

		for (Segment segment : segments) {
			PointKey p1 = new PointKey(segment.p1().getX(), segment.p1().getY());
			PointKey p2 = new PointKey(segment.p2().getX(), segment.p2().getY());
			skeletonEndpoints.add(p1);
			skeletonEndpoints.add(p2);
			skeletonLines.add(new SegmentLine(p1, p2));
		}

		for (PointKey v : inputVertices) {
			assertFalse(!skeletonEndpoints.contains(v), "Input vertex " + v + " not found in skeleton endpoints for " + csvFile.getFileName());
		}

		HPRtree tree = new HPRtree();
		for (int i = 0; i < skeletonLines.size(); i++) {
			tree.insert(skeletonLines.get(i).envelope(), Integer.valueOf(i));
		}
		tree.build();

		for (int i = 0; i < skeletonLines.size(); i++) {
			SegmentLine s1 = skeletonLines.get(i);
			@SuppressWarnings("unchecked")
			List<Integer> candidates = tree.query(s1.envelope());
			for (Integer candidate : candidates) {
				int j = candidate.intValue();
				if (j <= i) {
					continue;
				}
				SegmentLine s2 = skeletonLines.get(j);
				assertFalse(segmentsIntersectStrict(s1.p1, s1.p2, s2.p1, s2.p2),
						"Skeleton segments intersect in " + csvFile.getFileName() + ": " + s1 + " and " + s2);
			}
		}

		var faces = skeleton.asPolygonFaces();
		assertEquals(polygon.getArea(), faces.getArea(), 1e-9);
		var polys = (List<Polygon>) PolygonExtracter.getPolygons(faces);
		var polyArray = polys.toArray(new Geometry[0]);
		assertTrue(CoverageValidator.isValid(polyArray));
	}

	private static Stream<Path> csvFiles() {
		if (!Files.isDirectory(CSV_DIR)) {
			return Stream.empty();
		}
		try {
			return Files.list(CSV_DIR).filter(path -> path.getFileName().toString().toLowerCase().endsWith(".csv")).sorted().toList().stream();
		} catch (IOException e) {
			throw new RuntimeException("Failed to list CSV files under " + CSV_DIR, e);
		}
	}

	private static List<List<Coordinate>> readCsvPolygon(Path path) throws IOException {
		List<List<Coordinate>> rings = new ArrayList<>();
		for (String line : Files.readAllLines(path)) {
			if (line == null || line.isBlank()) {
				continue;
			}
			String[] tokens = line.split(",");
			List<Double> values = new ArrayList<>();
			for (String token : tokens) {
				String trimmed = token.trim();
				if (!trimmed.isEmpty()) {
					values.add(Double.parseDouble(trimmed));
				}
			}
			List<Coordinate> points = new ArrayList<>();
			for (int i = 0; i + 1 < values.size(); i += 2) {
				points.add(c(values.get(i), values.get(i + 1)));
			}
			if (!points.isEmpty()) {
				Coordinate first = points.get(0);
				Coordinate last = points.get(points.size() - 1);
				if (!first.equals2D(last)) {
					points.add(c(first.x, first.y));
				}
				rings.add(points);
			}
		}
		return rings;
	}

	private static Polygon toPolygon(List<List<Coordinate>> rings) {
		LinearRing shell = toLinearRing(rings.get(0));
		LinearRing[] holes = new LinearRing[Math.max(0, rings.size() - 1)];
		for (int i = 1; i < rings.size(); i++) {
			holes[i - 1] = toLinearRing(rings.get(i));
		}
		return GEOMETRY_FACTORY.createPolygon(shell, holes);
	}

	private static LinearRing toLinearRing(List<Coordinate> ringCoords) {
		return GEOMETRY_FACTORY.createLinearRing(ringCoords.toArray(Coordinate[]::new));
	}

	private static boolean segmentsIntersectStrict(PointKey p1, PointKey p2, PointKey p3, PointKey p4) {
		if (pointEqual(p1, p3) || pointEqual(p1, p4) || pointEqual(p2, p3) || pointEqual(p2, p4)) {
			return false;
		}

		int o1 = orientation(p1, p2, p3);
		int o2 = orientation(p1, p2, p4);
		int o3 = orientation(p3, p4, p1);
		int o4 = orientation(p3, p4, p2);

		return o1 != o2 && o3 != o4;
	}

	private static boolean pointEqual(PointKey a, PointKey b) {
		double tol = 1e-10;
		return Math.abs(a.x - b.x) < tol && Math.abs(a.y - b.y) < tol;
	}

	private static int orientation(PointKey p, PointKey q, PointKey r) {
		double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
		double tol = 1e-12;
		if (Math.abs(val) < tol) {
			return 0;
		}
		return val > 0 ? 1 : 2;
	}

	private static Coordinate c(double x, double y) {
		return new Coordinate(x, y);
	}

	private static Set<String> canonicalSegments(List<Segment> segments) {
		return segments.stream().map(s -> canonicalSegment(s.p1().getX(), s.p1().getY(), s.p2().getX(), s.p2().getY())).collect(Collectors.toSet());
	}

	private static String canonicalSegment(double x1, double y1, double x2, double y2) {
		String a = canonicalPoint(x1, y1);
		String b = canonicalPoint(x2, y2);
		return (a.compareTo(b) <= 0) ? (a + "|" + b) : (b + "|" + a);
	}

	private static String canonicalPoint(double x, double y) {
		return round6(x) + "," + round6(y);
	}

	private static double round6(double v) {
		return Math.rint(v * 1_000_000.0) / 1_000_000.0;
	}

	private static Integer findVertexInfo(InputMesh mesh, double x, double y) {
		for (var v : mesh.vertices) {
			if (pointEqual(new PointKey(v.x, v.y), new PointKey(x, y))) {
				return v.info;
			}
		}
		fail("Expected boundary vertex at " + x + "," + y);
		return null; // unreachable
	}

	private static Map<String, Integer> expectedSegments() {
		Map<String, Integer> expected = new HashMap<>();
		expected.put("eberly-10.csv", 36);
		expected.put("eberly-14.csv", 39);
		expected.put("elgindy-1.csv", 31);
		expected.put("gray-embroidery.csv", 61);
		expected.put("held-1.csv", 122);
		expected.put("held-12.csv", 63);
		expected.put("held-3.csv", 118);
		expected.put("held-7a.csv", 125);
		expected.put("held-7b.csv", 123);
		expected.put("held-7c.csv", 125);
		expected.put("held-7d.csv", 119);
		expected.put("mapbox-building.csv", 26);
		expected.put("mapbox-dude.csv", 211);
		expected.put("matisse-alga.csv", 494);
		expected.put("matisse-blue.csv", 394);
		expected.put("matisse-icarus.csv", 178);
		expected.put("matisse-nuit.csv", 245);
		expected.put("mei-2.csv", 49);
		expected.put("mei-3.csv", 81);
		expected.put("mei-4.csv", 231);
		expected.put("mei-5.csv", 555);
		expected.put("mei-6.csv", 732);
		expected.put("meisters-3.csv", 57);
		expected.put("misc-discobolus.csv", 291);
		expected.put("misc-fu.csv", 352);
		expected.put("seidel-3.csv", 37);
		expected.put("skimage-horse.csv", 235);
		expected.put("toussaint-1a.csv", 247);
		return expected;
	}

	private record PointKey(double x, double y) {
	}

	private record SegmentLine(PointKey p1, PointKey p2) {
		Envelope envelope() {
			return new Envelope(Math.min(p1.x, p2.x), Math.max(p1.x, p2.x), Math.min(p1.y, p2.y), Math.max(p1.y, p2.y));
		}
	}
}
