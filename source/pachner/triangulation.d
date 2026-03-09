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

import std.algorithm : canFind, max, reduce, remove, sort;
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

    /// f-vector: counts of distinct simplices of each dimension
    size_t nVertices;
    size_t nEdges;
    size_t nTriangles;
    size_t nTetrahedra;

    /// Add a tetrahedron with the given vertex labels
    size_t addTetrahedron(VertexLabel[4] verts)
    {
        size_t idx = tets.length;
        tets ~= Tet(verts);
        return idx;
    }

    /// Number of tetrahedra
    size_t size() const { return tets.length; }

    /// Recompute f-vector counts from scratch by enumerating all
    /// distinct vertices, edges, and triangles across all tetrahedra.
    void recount()
    {
        nTetrahedra = tets.length;

        // Count distinct vertices
        VertexLabel[] verts;
        foreach (ref t; tets)
            foreach (v; t.vertices)
                if (!canFind(verts, v))
                    verts ~= v;
        nVertices = verts.length;

        // Count distinct edges (sorted pairs)
        VertexLabel[2][] edgeList;
        foreach (ref t; tets)
            foreach (i; 0 .. 4)
                foreach (j; i + 1 .. 4)
                {
                    VertexLabel[2] e = [t.vertices[i], t.vertices[j]];
                    sort(e[]);
                    bool found = false;
                    foreach (ref existing; edgeList)
                        if (existing == e) { found = true; break; }
                    if (!found)
                        edgeList ~= e;
                }
        nEdges = edgeList.length;

        // Count distinct triangles (sorted triples)
        VertexLabel[3][] triList;
        foreach (ref t; tets)
            foreach (v; 0 .. 4)
            {
                auto f = t.face(cast(ubyte) v);
                bool found = false;
                foreach (ref existing; triList)
                    if (existing == f) { found = true; break; }
                if (!found)
                    triList ~= f;
            }
        nTriangles = triList.length;
    }

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

        // f-vector delta: +1 vertex, +4 edges, +6 triangles, +3 tets
        nVertices += 1;
        nEdges += 4;
        nTriangles += 6;
        nTetrahedra += 3;

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

        // f-vector delta: -1 vertex, -4 edges, -6 triangles, -3 tets
        nVertices -= 1;
        nEdges -= 4;
        nTriangles -= 6;
        nTetrahedra -= 3;

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

        // f-vector delta: +0 vertices, +1 edge, +2 triangles, +1 tet
        nEdges += 1;
        nTriangles += 2;
        nTetrahedra += 1;

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

        // f-vector delta: +0 vertices, -1 edge, -2 triangles, -1 tet
        nEdges -= 1;
        nTriangles -= 2;
        nTetrahedra -= 1;

        return true;
    }

    /// Returns a vertex label not currently used by any tetrahedron.
    /// Only available for integral VertexLabel types.
    static if (__traits(isIntegral, VertexLabel))
    VertexLabel nextVertexLabel() const
    {
        if (tets.length == 0)
            return 0;
        VertexLabel maxLabel = tets[0].vertices[0];
        foreach (ref t; tets)
            foreach (v; t.vertices)
                if (v > maxLabel)
                    maxLabel = v;
        return cast(VertexLabel)(maxLabel + 1);
    }

    /// Collect all distinct vertex labels used in the triangulation
    VertexLabel[] vertexLabels() const
    {
        VertexLabel[] labels;
        foreach (ref t; tets)
            foreach (v; t.vertices)
                if (!canFind(labels, v))
                    labels ~= v;
        sort(labels);
        return labels;
    }

    /// Collect all distinct edges (pairs of vertex labels) in the triangulation
    VertexLabel[2][] edges() const
    {
        VertexLabel[2][] result;
        foreach (ref t; tets)
        {
            foreach (i; 0 .. 4)
                foreach (j; i + 1 .. 4)
                {
                    VertexLabel[2] edge = [t.vertices[i], t.vertices[j]];
                    sort(edge[]);
                    bool found = false;
                    foreach (ref e; result)
                        if (e == edge) { found = true; break; }
                    if (!found)
                        result ~= edge;
                }
        }
        return result;
    }

    /**
     * Returns all valid Pachner moves for the current triangulation.
     *
     * Each move is represented as a PachnerMove tagged union that
     * records the move type and its parameters.
     */
    PachnerMove!VertexLabel[] validMoves() const
    {
        PachnerMove!VertexLabel[] moves;

        // 1-4 moves: one per tetrahedron (always valid since we use a fresh vertex)
        // Only available for integral vertex label types
        static if (__traits(isIntegral, VertexLabel))
        {
            VertexLabel fresh = nextVertexLabel();
            foreach (i; 0 .. tets.length)
            {
                moves ~= PachnerMove!VertexLabel.make14(i, fresh);
            }
        }

        // 2-3 moves: for each pair of tets sharing exactly 3 vertices
        foreach (i; 0 .. tets.length)
            foreach (j; i + 1 .. tets.length)
            {
                auto v1 = tets[i].vertices;
                auto v2 = tets[j].vertices;

                VertexLabel[] common;
                foreach (u; v1)
                    if (canFind(v2[], u) && !canFind(common, u))
                        common ~= u;

                if (common.length != 3)
                    continue;

                // Find non-shared vertices
                VertexLabel d, e;
                foreach (u; v1)
                    if (!canFind(common, u)) d = u;
                foreach (u; v2)
                    if (!canFind(common, u)) e = u;

                // Check validity: edge [d,e] must not exist in surviving tets
                if (!hasSimplexContaining([d, e], [i, j]))
                    moves ~= PachnerMove!VertexLabel.make23(i, j, common[0 .. 3], d, e);
            }

        // 3-2 moves: for each edge with exactly 3 tets around it
        auto allEdges = edges();
        foreach (ref edge; allEdges)
        {
            size_t[] tetIndices;
            foreach (i, ref t; tets)
                if (canFind(t.vertices[], edge[0]) && canFind(t.vertices[], edge[1]))
                    tetIndices ~= i;

            if (tetIndices.length != 3)
                continue;

            // Collect outer vertices
            VertexLabel[] outer;
            foreach (idx; tetIndices)
                foreach (u; tets[idx].vertices)
                    if (u != edge[0] && u != edge[1] && !canFind(outer, u))
                        outer ~= u;

            if (outer.length != 3)
                continue;

            // Verify each tet has exactly 2 outer vertices
            bool valid = true;
            foreach (idx; tetIndices)
            {
                int outerCount = 0;
                foreach (u; tets[idx].vertices)
                    if (canFind(outer, u))
                        outerCount++;
                if (outerCount != 2) { valid = false; break; }
            }
            if (!valid) continue;

            // Check validity: face [a,b,c] must not exist in surviving tets
            if (!hasSimplexContaining(outer, tetIndices))
                moves ~= PachnerMove!VertexLabel.make32(edge, outer[0 .. 3]);
        }

        // 4-1 moves: for each vertex with exactly 4 tets around it
        auto allVerts = vertexLabels();
        foreach (v; allVerts)
        {
            size_t[] tetIndices;
            foreach (i, ref t; tets)
                if (canFind(t.vertices[], v))
                    tetIndices ~= i;

            if (tetIndices.length != 4)
                continue;

            VertexLabel[] outer;
            foreach (idx; tetIndices)
                foreach (u; tets[idx].vertices)
                    if (u != v && !canFind(outer, u))
                        outer ~= u;

            if (outer.length != 4)
                continue;

            bool valid = true;
            foreach (idx; tetIndices)
            {
                int outerCount = 0;
                foreach (u; tets[idx].vertices)
                    if (canFind(outer, u))
                        outerCount++;
                if (outerCount != 3) { valid = false; break; }
            }
            if (!valid) continue;

            moves ~= PachnerMove!VertexLabel.make41(v, outer[0 .. 4]);
        }

        return moves;
    }

    /// Return a deep copy of this triangulation
    Triangulation dup() const
    {
        Triangulation copy;
        copy.tets = tets.dup;
        copy.nVertices = nVertices;
        copy.nEdges = nEdges;
        copy.nTriangles = nTriangles;
        copy.nTetrahedra = nTetrahedra;
        return copy;
    }

    /// Pretty-print summary
    void print() const
    {
        writefln("Triangulation: %d vertices, %d edges, %d triangles, %d tetrahedra",
            nVertices, nEdges, nTriangles, nTetrahedra);
        foreach (i, ref t; tets)
        {
            writefln("  Tet %d: %s", i, t.vertices);
        }
    }
}

/**
 * Tagged union representing a Pachner move and its parameters.
 * Stores enough information to perform the move and compute its inverse.
 */
struct PachnerMove(VertexLabel = size_t)
{
    enum Type { move14, move41, move23, move32 }
    Type type;

    // 1-4 / 4-1 parameters
    size_t tetIdx;        // tet index for 1-4
    VertexLabel newVertex; // new vertex for 1-4, vertex to remove for 4-1
    VertexLabel[4] outerVertices; // outer vertices for 4-1

    // 2-3 / 3-2 parameters
    size_t tetIdx1, tetIdx2;     // tet indices for 2-3
    VertexLabel[3] faceVertices; // shared face for 2-3, outer vertices for 3-2
    VertexLabel[2] edgeVertices; // edge [d,e] for 2-3, shared edge for 3-2

    static PachnerMove make14(size_t tetIdx, VertexLabel newVertex)
    {
        PachnerMove m;
        m.type = Type.move14;
        m.tetIdx = tetIdx;
        m.newVertex = newVertex;
        return m;
    }

    static PachnerMove make41(VertexLabel v, VertexLabel[4] outer)
    {
        PachnerMove m;
        m.type = Type.move41;
        m.newVertex = v;
        m.outerVertices = outer;
        return m;
    }

    static PachnerMove make23(size_t i, size_t j, VertexLabel[3] face, VertexLabel d, VertexLabel e)
    {
        PachnerMove m;
        m.type = Type.move23;
        m.tetIdx1 = i;
        m.tetIdx2 = j;
        m.faceVertices = face;
        m.edgeVertices = [d, e];
        return m;
    }

    static PachnerMove make32(VertexLabel[2] edge, VertexLabel[3] outer)
    {
        PachnerMove m;
        m.type = Type.move32;
        m.edgeVertices = edge;
        m.faceVertices = outer;
        return m;
    }

    /// Returns the inverse move that undoes this one.
    PachnerMove inverse() const
    {
        PachnerMove inv;
        final switch (type)
        {
            case Type.move14:
                // Inverse of 1-4 is 4-1 on the new vertex
                inv.type = Type.move41;
                inv.newVertex = newVertex;
                // outerVertices not known here; will be resolved at execution
                break;
            case Type.move41:
                // Inverse of 4-1 is 1-4 on the resulting tet
                inv.type = Type.move14;
                inv.newVertex = newVertex;
                inv.outerVertices = outerVertices;
                // tetIdx will need to be found at execution
                break;
            case Type.move23:
                // Inverse of 2-3 is 3-2 on edge [d,e]
                inv.type = Type.move32;
                inv.edgeVertices = edgeVertices;
                inv.faceVertices = faceVertices;
                break;
            case Type.move32:
                // Inverse of 3-2 is 2-3 on the two tets sharing face
                inv.type = Type.move23;
                inv.edgeVertices = edgeVertices;
                inv.faceVertices = faceVertices;
                // tetIdx1/2 will need to be found at execution
                break;
        }
        return inv;
    }

    /// Find the index of a tet matching the given sorted vertices, or size_t.max
    private static size_t findTet(ref const Triangulation!VertexLabel tri, VertexLabel[4] sorted)
    {
        foreach (i, ref t; tri.tets)
            if (t.sortedVertices() == sorted)
                return i;
        return size_t.max;
    }

    /// Execute this move on the given triangulation. Returns true on success.
    bool execute(ref Triangulation!VertexLabel tri)
    {
        final switch (type)
        {
            case Type.move14:
                // Find the tet by vertex content if outerVertices is set
                size_t idx = tetIdx;
                if (outerVertices != typeof(outerVertices).init)
                {
                    VertexLabel[4] sorted = outerVertices;
                    sort(sorted[]);
                    idx = findTet(tri, sorted);
                    if (idx == size_t.max) return false;
                }
                auto result = tri.move14(idx, newVertex);
                return result[0] != size_t.max;

            case Type.move41:
                return tri.move41(newVertex);

            case Type.move23:
                // Find the two tets by their vertex content
                VertexLabel[4] tet1verts, tet2verts;
                foreach (i; 0 .. 3) tet1verts[i] = faceVertices[i];
                tet1verts[3] = edgeVertices[0];
                foreach (i; 0 .. 3) tet2verts[i] = faceVertices[i];
                tet2verts[3] = edgeVertices[1];
                sort(tet1verts[]);
                sort(tet2verts[]);

                auto idx1 = findTet(tri, tet1verts);
                auto idx2 = findTet(tri, tet2verts);
                if (idx1 == size_t.max || idx2 == size_t.max)
                    return false;
                return tri.move23(idx1, idx2);

            case Type.move32:
                return tri.move32(edgeVertices[0], edgeVertices[1]);
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

    auto original = tri.dup;
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

    auto original = tri.dup;

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

    auto original = tri.dup;

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

/// Build the boundary of the 4-simplex [0,1,2,3,4] — a triangulation of S^3
Triangulation!size_t fourSimplexBoundary()
{
    Triangulation!size_t tri;
    // 5 tetrahedra, each omitting one vertex from {0,1,2,3,4}
    tri.addTetrahedron([1, 2, 3, 4]);
    tri.addTetrahedron([0, 2, 3, 4]);
    tri.addTetrahedron([0, 1, 3, 4]);
    tri.addTetrahedron([0, 1, 2, 4]);
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.recount();
    return tri;
}

unittest
{
    // Verify 4-simplex boundary has expected structure
    auto tri = fourSimplexBoundary();
    assert(tri.size == 5);
    assert(tri.vertexLabels() == [0, 1, 2, 3, 4]);

    // f-vector of boundary of 4-simplex: 5 vertices, 10 edges, 10 triangles, 5 tets
    assert(tri.nVertices == 5);
    assert(tri.nEdges == 10);
    assert(tri.nTriangles == 10);
    assert(tri.nTetrahedra == 5);

    // Should have valid moves available
    auto moves = tri.validMoves();
    assert(moves.length > 0);
}

unittest
{
    // f-vector is maintained correctly through 1-4 and 4-1 moves
    auto tri = fourSimplexBoundary();

    tri.move14(0, 5);
    assert(tri.nVertices == 6);
    assert(tri.nEdges == 14);
    assert(tri.nTriangles == 16);
    assert(tri.nTetrahedra == 8);
    tri.recount();
    assert(tri.nVertices == 6);
    assert(tri.nEdges == 14);
    assert(tri.nTriangles == 16);
    assert(tri.nTetrahedra == 8);

    tri.move41(5);
    assert(tri.nVertices == 5);
    assert(tri.nEdges == 10);
    assert(tri.nTriangles == 10);
    assert(tri.nTetrahedra == 5);
}

unittest
{
    // f-vector is maintained correctly through 2-3 and 3-2 moves
    Triangulation!size_t tri;
    tri.addTetrahedron([0, 1, 2, 3]);
    tri.addTetrahedron([0, 1, 2, 4]);
    tri.recount();
    assert(tri.nVertices == 5);
    assert(tri.nEdges == 9);
    assert(tri.nTriangles == 7);
    assert(tri.nTetrahedra == 2);

    tri.move23(0, 1);
    assert(tri.nVertices == 5);
    assert(tri.nEdges == 10);
    assert(tri.nTriangles == 9);
    assert(tri.nTetrahedra == 3);
    // Verify against recount
    tri.recount();
    assert(tri.nEdges == 10);
    assert(tri.nTriangles == 9);

    tri.move32(3, 4);
    assert(tri.nVertices == 5);
    assert(tri.nEdges == 9);
    assert(tri.nTriangles == 7);
    assert(tri.nTetrahedra == 2);
}

unittest
{
    import std.random : Random, uniform;

    // Verify f-vector stays consistent with recount() through random walk
    auto tri = fourSimplexBoundary();
    auto rng = Random(123);
    alias Move = PachnerMove!size_t;

    foreach (step; 0 .. 50)
    {
        auto moves = tri.validMoves();
        auto idx = uniform(0, moves.length, rng);
        auto move = moves[idx];
        if (move.type == Move.Type.move14)
            move.newVertex = tri.nextVertexLabel();
        move.execute(tri);

        // Verify tracked counts match full recount
        auto v = tri.nVertices;
        auto e = tri.nEdges;
        auto f = tri.nTriangles;
        auto t = tri.nTetrahedra;
        tri.recount();
        assert(tri.nVertices == v);
        assert(tri.nEdges == e);
        assert(tri.nTriangles == f);
        assert(tri.nTetrahedra == t);
    }
}

unittest
{
    import std.random : Random, uniform;

    // Random walk of 100 Pachner moves starting from the 4-simplex boundary,
    // then undo all moves in reverse to recover the original.
    auto tri = fourSimplexBoundary();
    auto original = tri.dup;

    alias Move = PachnerMove!size_t;
    Move[] history;
    auto rng = Random(42); // fixed seed for reproducibility

    // Perform 100 random moves
    foreach (step; 0 .. 100)
    {
        auto moves = tri.validMoves();
        assert(moves.length > 0, "No valid moves available at step " ~ (cast(int) step).stringof);

        // Pick a random valid move
        auto idx = uniform(0, moves.length, rng);
        auto move = moves[idx];

        // For 1-4 moves, use nextVertexLabel for the new vertex
        if (move.type == Move.Type.move14)
            move.newVertex = tri.nextVertexLabel();

        bool ok = move.execute(tri);
        assert(ok, "Move execution failed at step " ~ (cast(int) step).stringof);
        history ~= move;
    }

    // Undo all moves in reverse order
    foreach_reverse (ref move; history)
    {
        auto inv = move.inverse();
        bool ok = inv.execute(tri);
        assert(ok, "Inverse move failed");
    }

    assert(tri == original, "Triangulation not recovered after undoing all moves");
}
