# Clipper2 test data

These files are copied unchanged from Clipper2 2.0.1 (`Tests/` in
https://github.com/AngusJohnson/Clipper2), under the Boost Software License
in `LICENSE` here.

Test 85 (`tests/85_clipper2_parity.c`) runs `Polygons.txt` and `Lines.txt`
through `qaws_clip` with the tolerances of Clipper2's own tests. It also
checks every result against the definition, on a grid of points per record.
