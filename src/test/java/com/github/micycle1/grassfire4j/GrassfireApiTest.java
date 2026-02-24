package com.github.micycle1.grassfire4j;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.DynamicTest.dynamicTest;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestFactory;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.strtree.STRtree;

import com.github.micycle1.grassfire4j.model.Model.Skeleton.Segment;

class GrassfireApiTest {

	private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();
	private static final Path CSV_DIR = Path.of("src", "test", "resources", "csv");
	private static final Map<String, Integer> EXPECTED_SEGMENTS = expectedSegments();

	@Test
	void internalSegmentsCount() {
		List<List<Coordinate>> rings = List.of(List.of(
				c(0.0, 0.0),
				c(20.0, 0.0),
				c(20.0, 10.0),
				c(10.0, 10.0),
				c(10.0, 20.0),
				c(0.0, 20.0),
				c(0.0, 0.0)));

		Polygon polygon = toPolygon(rings);
		var skeleton = GrassfireApi.computeSkeleton(polygon);
		assertEquals(8, skeleton.segments().size(), "Unexpected internal segment count");
	}

	@TestFactory
	Stream<DynamicTest> skeletonIntegrity() throws IOException {
		return csvFiles().map(csvFile -> dynamicTest(csvFile.getFileName().toString(), () -> runSkeletonIntegrity(csvFile)));
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
		var skeleton = GrassfireApi.computeSkeleton(polygon);
		List<Segment> segments = skeleton.segments();

		Integer expected = EXPECTED_SEGMENTS.get(csvFile.getFileName().toString());
		if (expected != null) {
			assertEquals(expected.intValue(), segments.size(),
					"Unexpected segment count for " + csvFile.getFileName());
		}

		Set<PointKey> skeletonEndpoints = new HashSet<>();
		List<SegmentLine> skeletonLines = new ArrayList<>(segments.size());

		for (Segment segment : segments) {
			PointKey p1 = new PointKey(segment.p1().x(), segment.p1().y());
			PointKey p2 = new PointKey(segment.p2().x(), segment.p2().y());
			skeletonEndpoints.add(p1);
			skeletonEndpoints.add(p2);
			skeletonLines.add(new SegmentLine(p1, p2));
		}

		for (PointKey v : inputVertices) {
			assertFalse(!skeletonEndpoints.contains(v),
					"Input vertex " + v + " not found in skeleton endpoints for " + csvFile.getFileName());
		}

		STRtree tree = new STRtree();
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
	}

	private static Stream<Path> csvFiles() {
		if (!Files.isDirectory(CSV_DIR)) {
			return Stream.empty();
		}
		try {
			return Files.list(CSV_DIR)
					.filter(path -> path.getFileName().toString().toLowerCase().endsWith(".csv"))
					.sorted()
					.toList()
					.stream();
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
		double tol = 1e-9;
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
			return new Envelope(
					Math.min(p1.x, p2.x),
					Math.max(p1.x, p2.x),
					Math.min(p1.y, p2.y),
					Math.max(p1.y, p2.y));
		}
	}
}
