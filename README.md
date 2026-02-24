[![](https://jitpack.io/v/micycle1/grassfire4j.svg)](https://jitpack.io/#micycle1/grassfire4j)

# grassfire4j

A kinetic straight‑skeleton (wavefront‑collapse / “grassfire”) implementation in Java.

## Overview

grassfire4j implements a kinetic straight‑skeleton algorithm. The implementation is rooted in the original Python [reference](https://github.com/bmmeijers/grassfire), which I rearchitected and improved in my redesigned Python [implementation](https://github.com/micycle1/grassfire2). That improved Python design – which itself incorporated useful ideas from the academic [surfer2](https://github.com/cgalab/surfer2) C++ project – was then ported and adapted to Java to make this project.

This project is the first full kinetic straight‑skeleton implementation available in Java. Compared to the other notable Java implementation, [campskeleton](https://github.com/twak/campskeleton), which uses Felkel’s edge‑collision approach, grassfire4j uses a triangulation + kinetic event queue approach with careful tie‑breaking and degenerate‑case handling. In practice this yields better behaviour and significantly improved performance on polygonal inputs.

## Features
- Kinetic (wavefront‑collapse) straight‑skeleton computed entirely in Java.
- Accepts JTS `Polygon` inputs; supports polygons with holes.
- Adapter-based input pipeline (`InputMesh`) makes it easy to plug in user-supplied adapters for other input types.
- Supports variable edge weights.
- Produces a Skeleton model with nodes, kinetic vertices and skeleton segments suitable for visualisation or export.
