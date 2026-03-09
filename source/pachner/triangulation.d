/**
 * Triangulation data structure for 3-manifolds.
 *
 * Stores a collection of tetrahedra whose gluings are determined
 * implicitly by shared vertex labels — two tetrahedra sharing 3
 * vertex labels are glued along that common face.
 *
 * Supports Pachner moves (bistellar flips):
 *   - 2-3 move: replace 2 tetrahedra sharing a face with 3 sharing an edge
 *   - 3-2 move: inverse of the above
 *   - 1-4 move: subdivide 1 tetrahedron into 4
 *   - 4-1 move: inverse of the above
 */
module pachner.triangulation;

import pachner.simplex;

import std.algorithm : sort;
import std.stdio : writefln;

/**
 * A triangulation of a 3-manifold stored as a list of tetrahedra.
 *
 * Gluings are implicit: tetrahedra that share 3 vertex labels
 * are considered glued along the corresponding face.
 */
struct Triangulation(VertexLabel = size_t)
{
    alias Tet = Tetrahedron!VertexLabel;

    /// All tetrahedra in the triangulation
    Tet[] tets;

    /// Add a tetrahedron with the given vertex labels
    size_t addTetrahedron(VertexLabel[4] verts)
    {
        size_t idx = tets.length;
        tets ~= Tet(verts);
        return idx;
    }

    /// Number of tetrahedra
    size_t size() const { return tets.length; }

    /// Pretty-print summary
    void print() const
    {
        writefln("Triangulation: %d tetrahedra", tets.length);
        foreach (i, ref t; tets)
        {
            writefln("  Tet %d: %s", i, t.vertices);
        }
    }
}

unittest
{
    Triangulation!size_t tri;

    tri.addTetrahedron([0, 1, 2, 3]);
    tri.addTetrahedron([0, 1, 2, 4]);
    assert(tri.size == 2);

    // The two tets share face [0, 1, 2] implicitly via shared vertex labels
}

unittest
{
    // Works with string labels too
    Triangulation!string tri;
    tri.addTetrahedron(["a", "b", "c", "d"]);
    assert(tri.size == 1);
}
