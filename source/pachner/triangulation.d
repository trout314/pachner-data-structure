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
import std.array : array;
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

    /// Order-independent equality: same set of tetrahedra (ignoring
    /// vertex order within each tet and tet order in the list)
    bool opEquals(ref const Triangulation rhs) const
    {
        if (tets.length != rhs.tets.length)
            return false;

        // Sort copies of both tet lists by sorted vertex labels
        auto a = tets.dup;
        auto b = rhs.tets.dup;
        sort(a);
        sort(b);

        foreach (i; 0 .. a.length)
            if (a[i].sortedVertices() != b[i].sortedVertices())
                return false;

        return true;
    }

    /// Returns true if any tet (not in excludeIndices) contains all
    /// the given vertex labels.
    bool hasSimplexContaining(const VertexLabel[] vertices, const size_t[] excludeIndices = null) const
    {
        foreach (i, ref t; tets)
        {
            if (excludeIndices !is null && canFind(excludeIndices, i))
                continue;
            bool allFound = true;
            foreach (v; vertices)
                if (!canFind(t.vertices[], v))
                    { allFound = false; break; }
            if (allFound)
                return true;
        }
        return false;
    }

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
        // Check that v is not already a vertex of any existing tet
        if (hasSimplexContaining([v]))
            return [size_t.max, size_t.max, size_t.max, size_t.max];

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

    /**
     * Perform a 2-3 Pachner move on two tetrahedra that share a face.
     *
     * Given tets [a,b,c,d] and [a,b,c,e] sharing face [a,b,c],
     * replaces them with 3 tets sharing edge (d,e):
     *   [a,b,d,e], [a,c,d,e], [b,c,d,e]
     *
     * The caller specifies the two tet indices. The move fails if
     * they don't share exactly 3 vertex labels.
     *
     * Returns: true if the move was performed, false if preconditions fail.
     */
    bool move23(size_t tetIdx1, size_t tetIdx2)
    {
        if (tetIdx1 >= tets.length || tetIdx2 >= tets.length)
            return false;
        if (tetIdx1 == tetIdx2)
            return false;

        auto v1 = tets[tetIdx1].vertices;
        auto v2 = tets[tetIdx2].vertices;

        // Find shared and non-shared vertices
        VertexLabel[] common;
        foreach (u; v1)
            if (canFind(v2[], u) && !canFind(common, u))
                common ~= u;

        if (common.length != 3)
            return false;

        // Find the non-shared vertex from each tet
        VertexLabel d, e;
        bool foundD = false, foundE = false;
        foreach (u; v1)
            if (!canFind(common, u)) { d = u; foundD = true; }
        foreach (u; v2)
            if (!canFind(common, u)) { e = u; foundE = true; }

        if (!foundD || !foundE)
            return false;

        // Check that no surviving tet contains edge [d, e]
        if (hasSimplexContaining([d, e], [tetIdx1, tetIdx2]))
            return false;

        // Remove both tets (higher index first)
        if (tetIdx1 > tetIdx2)
        {
            removeTetrahedron(tetIdx1);
            removeTetrahedron(tetIdx2);
        }
        else
        {
            removeTetrahedron(tetIdx2);
            removeTetrahedron(tetIdx1);
        }

        // Add 3 new tets, each with edge (d,e) and one edge of the shared triangle
        foreach (i; 0 .. 3)
        {
            VertexLabel[4] verts;
            size_t k = 0;
            foreach (j; 0 .. 3)
                if (j != i)
                    verts[k++] = common[j];
            verts[2] = d;
            verts[3] = e;
            addTetrahedron(verts);
        }

        return true;
    }

    /**
     * Perform a 3-2 Pachner move on 3 tetrahedra sharing edge (d, e).
     *
     * Finds all tets containing both d and e. If there are exactly 3,
     * and the 3 "other" vertices (excluding d and e) are all distinct,
     * replaces them with 2 tets:
     *   [a, b, c, d] and [a, b, c, e]
     * where {a, b, c} are the 3 outer vertices.
     *
     * Returns: true if the move was performed, false if preconditions fail.
     */
    bool move32(VertexLabel d, VertexLabel e)
    {
        // Find all tets containing both d and e
        size_t[] tetIndices;
        foreach (i, ref t; tets)
            if (canFind(t.vertices[], d) && canFind(t.vertices[], e))
                tetIndices ~= i;

        if (tetIndices.length != 3)
            return false;

        // Collect distinct non-{d,e} vertices
        VertexLabel[] outer;
        foreach (idx; tetIndices)
            foreach (u; tets[idx].vertices)
                if (u != d && u != e && !canFind(outer, u))
                    outer ~= u;

        if (outer.length != 3)
            return false;

        // Verify each tet has exactly 2 of the 3 outer vertices (plus d and e)
        foreach (idx; tetIndices)
        {
            int outerCount = 0;
            foreach (u; tets[idx].vertices)
                if (canFind(outer, u))
                    outerCount++;
            if (outerCount != 2)
                return false;
        }

        // Check that no surviving tet contains face [a, b, c]
        if (hasSimplexContaining(outer, tetIndices))
            return false;

        // Remove the 3 tets (highest index first)
        sort!"a > b"(tetIndices);
        foreach (idx; tetIndices)
            removeTetrahedron(idx);

        // Add 2 new tets
        addTetrahedron([outer[0], outer[1], outer[2], d]);
        addTetrahedron([outer[0], outer[1], outer[2], e]);

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
    // 4-1 move: round-trip with 1-4 recovers original triangulation
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);

    Triangulation!size_t original = tri;
    tri.move14(0, 4);
    assert(tri.size == 4);

    bool ok = tri.move41(4);
    assert(ok);
    assert(tri == original);
}

unittest
{
    // 4-1 move should fail if vertex doesn't have exactly 4 tets
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    assert(!tri.move41(0)); // vertex 0 is in only 1 tet
}

unittest
{
    // 2-3 move: two tets sharing a face become three
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.addTetrahedron([0, 1, 2, 4]);
    assert(tri.size == 2);

    bool ok = tri.move23(0, 1);
    assert(ok);
    assert(tri.size == 3);

    // Each new tet should contain both non-shared vertices (3 and 4)
    foreach (ref t; tri.tets)
    {
        assert(canFind(t.vertices[], cast(size_t) 3));
        assert(canFind(t.vertices[], cast(size_t) 4));
    }
}

unittest
{
    // 2-3 move should fail if tets don't share exactly 3 vertices
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.addTetrahedron([4, 5, 6, 7]);
    assert(!tri.move23(0, 1));
}

unittest
{
    // 3-2 move: three tets sharing an edge become two
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 3, 4]);
    tri.addTetrahedron([0, 2, 3, 4]);
    tri.addTetrahedron([1, 2, 3, 4]);
    assert(tri.size == 3);

    bool ok = tri.move32(3, 4);
    assert(ok);
    assert(tri.size == 2);

    // Each new tet should contain all 3 outer vertices (0, 1, 2)
    foreach (ref t; tri.tets)
    {
        assert(canFind(t.vertices[], cast(size_t) 0));
        assert(canFind(t.vertices[], cast(size_t) 1));
        assert(canFind(t.vertices[], cast(size_t) 2));
    }
}

unittest
{
    // 2-3 then 3-2 round-trip recovers original
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.addTetrahedron([0, 1, 2, 4]);

    Triangulation!size_t original = tri;

    tri.move23(0, 1);
    assert(tri.size == 3);

    bool ok = tri.move32(3, 4);
    assert(ok);
    assert(tri == original);
}

unittest
{
    // 3-2 then 2-3 round-trip recovers original
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 3, 4]);
    tri.addTetrahedron([0, 2, 3, 4]);
    tri.addTetrahedron([1, 2, 3, 4]);

    Triangulation!size_t original = tri;

    tri.move32(3, 4);
    assert(tri.size == 2);

    bool ok = tri.move23(0, 1);
    assert(ok);
    assert(tri == original);
}

unittest
{
    // 3-2 move should fail if edge doesn't have exactly 3 tets
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    assert(!tri.move32(2, 3)); // edge (2,3) is in only 1 tet
}

unittest
{
    // 1-4 move should fail if vertex already exists
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    auto result = tri.move14(0, 3); // vertex 3 already in use
    assert(result == [size_t.max, size_t.max, size_t.max, size_t.max]);
    assert(tri.size == 1); // triangulation unchanged
}

unittest
{
    // 2-3 move should fail if edge [d,e] already exists in a surviving tet
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.addTetrahedron([0, 1, 2, 4]);
    tri.addTetrahedron([0, 3, 4, 5]); // surviving tet contains edge [3,4]
    assert(!tri.move23(0, 1)); // would create edge [3,4] which already exists
    assert(tri.size == 3); // triangulation unchanged
}

unittest
{
    // 3-2 move should fail if face [a,b,c] already exists in a surviving tet
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 3, 4]);
    tri.addTetrahedron([0, 2, 3, 4]);
    tri.addTetrahedron([1, 2, 3, 4]);
    tri.addTetrahedron([0, 1, 2, 5]); // surviving tet contains face [0,1,2]
    assert(!tri.move32(3, 4)); // would create face [0,1,2] which already exists
    assert(tri.size == 4); // triangulation unchanged
}

unittest
{
    // opEquals: same tets in different order are equal
    Triangulation!size_t a, b;
    a.addTetrahedron([0, 1, 2, 3]);
    a.addTetrahedron([0, 1, 2, 4]);

    b.addTetrahedron([0, 1, 2, 4]);
    b.addTetrahedron([0, 1, 2, 3]);
    assert(a == b);
}

unittest
{
    // opEquals: same tet with different vertex order are equal
    Triangulation!size_t a, b;
    a.addTetrahedron([0, 1, 2, 3]);
    b.addTetrahedron([3, 2, 1, 0]);
    assert(a == b);
}

unittest
{
    // opEquals: different tets are not equal
    Triangulation!size_t a, b;
    a.addTetrahedron([0, 1, 2, 3]);
    b.addTetrahedron([0, 1, 2, 4]);
    assert(a != b);
}
