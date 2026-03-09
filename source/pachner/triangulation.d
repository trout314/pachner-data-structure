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

import std.algorithm : canFind, remove, sort;
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

    /**
     * Perform a 4-1 Pachner move, removing vertex v.
     *
     * Finds all tetrahedra containing v. If there are exactly 4 and
     * they share exactly 4 distinct non-v vertices (each tet using 3
     * of them), removes the 4 tetrahedra and replaces them with a
     * single tetrahedron on the 4 outer vertices.
     *
     * Returns: true if the move was performed, false if preconditions fail.
     */
    bool move41(VertexLabel v)
    {
        // Find all tetrahedra containing v
        size_t[] tetIndices;
        foreach (i, ref t; tets)
            if (canFind(t.vertices[], v))
                tetIndices ~= i;

        if (tetIndices.length != 4)
            return false;

        // Collect all distinct non-v vertices
        VertexLabel[] outer;
        foreach (idx; tetIndices)
            foreach (u; tets[idx].vertices)
                if (u != v && !canFind(outer, u))
                    outer ~= u;

        if (outer.length != 4)
            return false;

        // Verify each tet has exactly 3 of the outer vertices (plus v)
        foreach (idx; tetIndices)
        {
            int outerCount = 0;
            foreach (u; tets[idx].vertices)
                if (canFind(outer, u))
                    outerCount++;
            if (outerCount != 3)
                return false;
        }

        // Remove the 4 tets (remove from highest index first to preserve indices)
        sort!"a > b"(tetIndices);
        foreach (idx; tetIndices)
            removeTetrahedron(idx);

        // Add the single replacement tetrahedron
        addTetrahedron(outer[0 .. 4]);

        return true;
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

unittest
{
    // 4-1 move: four tets become one (round-trip with 1-4)
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.move14(0, 4);
    assert(tri.size == 4);

    bool ok = tri.move41(4);
    assert(ok);
    assert(tri.size == 1);

    // The remaining tet should have exactly the original 4 vertices
    VertexLabel[] verts;
    foreach (u; tri.tets[0].vertices)
        if (!canFind(verts, u))
            verts ~= u;
    sort(verts);
    assert(verts == [0, 1, 2, 3]);
}

unittest
{
    // 4-1 move should fail if vertex doesn't have exactly 4 tets
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    assert(!tri.move41(0)); // vertex 0 is in only 1 tet
}

private alias VertexLabel = size_t; // for unittest visibility
