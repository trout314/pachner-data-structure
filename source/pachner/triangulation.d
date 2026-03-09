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

import std.algorithm : remove, sort;
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

    /// Remove tetrahedron at the given index (swaps with last element)
    void removeTetrahedron(size_t idx)
    in (idx < tets.length)
    {
        tets[idx] = tets[$ - 1];
        tets = tets[0 .. $ - 1];
    }

    /**
     * Perform a 1-4 Pachner move on the tetrahedron at index tetIdx.
     *
     * Replaces the tetrahedron [a, b, c, d] with 4 new tetrahedra,
     * each formed by replacing one original vertex with the new vertex v:
     *   [v, b, c, d], [a, v, c, d], [a, b, v, d], [a, b, c, v]
     *
     * The caller must provide a new vertex label v that is not already
     * used in the triangulation.
     *
     * Returns: the indices of the 4 new tetrahedra.
     */
    size_t[4] move14(size_t tetIdx, VertexLabel v)
    in (tetIdx < tets.length)
    {
        VertexLabel[4] orig = tets[tetIdx].vertices;
        removeTetrahedron(tetIdx);

        size_t[4] newIndices;
        foreach (i; 0 .. 4)
        {
            VertexLabel[4] verts = orig;
            verts[i] = v;
            newIndices[i] = addTetrahedron(verts);
        }
        return newIndices;
    }

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

unittest
{
    // 1-4 move: one tet becomes four
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    assert(tri.size == 1);

    auto newIdx = tri.move14(0, 4);
    assert(tri.size == 4);

    // Each new tet should contain vertex 4
    foreach (i; newIdx)
    {
        bool hasNew = false;
        foreach (v; tri.tets[i].vertices)
            if (v == 4) hasNew = true;
        assert(hasNew);
    }

    // Each new tet should have exactly 3 of the original vertices
    foreach (i; newIdx)
    {
        int origCount = 0;
        foreach (v; tri.tets[i].vertices)
            foreach (orig; [0, 1, 2, 3])
                if (v == orig) origCount++;
        assert(origCount == 3);
    }
}
