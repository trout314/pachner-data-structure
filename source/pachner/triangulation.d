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

import std.algorithm : canFind, countUntil, max, reduce, remove, sort, SwapStrategy;
import std.array : array;
import std.stdio : writefln;
import std.typecons : Nullable, nullable;

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

    /// Degree maps: number of tetrahedra containing each vertex/edge
    size_t[VertexLabel] vertexDegrees;
    size_t[VertexLabel[2]] edgeDegrees;

    /// Face index: face (sorted triple) → list of tet indices containing that face
    size_t[][VertexLabel[3]] faceTets;

    /// Parallel arrays for O(1) random access into index keys
    VertexLabel[] vertexKeys;
    size_t[VertexLabel] vertexKeyIndex;
    VertexLabel[2][] edgeKeys;
    size_t[VertexLabel[2]] edgeKeyIndex;
    VertexLabel[3][] faceKeys;
    size_t[VertexLabel[3]] faceKeyIndex;

    /// Running sums for variance computation (avoid full iteration)
    long vertexDegreeSum;
    long vertexDegreeSqSum;
    long edgeDegreeSum;
    long edgeDegreeSqSum;

    /// Monotonically increasing counter for fresh vertex labels
    static if (__traits(isIntegral, VertexLabel))
        VertexLabel nextVertexCounter = 0;

    /// Add a key to a parallel key array with swap-with-last support
    private static void addKey(K)(ref K[] arr, ref size_t[K] index, K key)
    {
        index[key] = arr.length;
        arr ~= key;
    }

    /// Remove a key from a parallel key array using swap-with-last
    private static void removeKey(K)(ref K[] arr, ref size_t[K] index, K key)
    {
        auto pos = index[key];
        auto lastPos = arr.length - 1;
        if (pos != lastPos)
        {
            auto last = arr[lastPos];
            arr[pos] = last;
            index[last] = pos;
        }
        arr = arr[0 .. lastPos];
        index.remove(key);
    }

    /// Add a tetrahedron with the given vertex labels.
    /// Updates all indexes (face, degree maps, key arrays).
    size_t addTetrahedron(VertexLabel[4] verts)
    {
        size_t idx = tets.length;
        tets ~= Tet(verts);

        // Update face index
        auto tet = Tet(verts);
        foreach (fi; 0 .. 4)
        {
            auto fk = tet.face(cast(ubyte) fi);
            if (fk !in faceTets)
                addKey(faceKeys, faceKeyIndex, fk);
            faceTets[fk] ~= idx;
        }

        // Update degree maps
        adjustDegrees(verts, +1);
        nTriangles = faceTets.length;
        nTetrahedra = tets.length;

        return idx;
    }

    /// Number of tetrahedra
    size_t size() const { return tets.length; }

    /// Return a sorted edge key for use in the edgeDegrees AA
    static VertexLabel[2] edgeKey(VertexLabel a, VertexLabel b)
    {
        VertexLabel[2] e = [a, b];
        sort(e[]);
        return e;
    }

    /// Update degree tracking when a tetrahedron is added (+1) or removed (-1)
    private void adjustDegrees(VertexLabel[4] verts, int delta)
    {
        // Update vertex degrees
        foreach (v; verts)
        {
            auto p = v in vertexDegrees;
            long oldDeg = p ? cast(long) *p : 0;
            long newDeg = oldDeg + delta;

            vertexDegreeSum += delta;
            vertexDegreeSqSum += newDeg * newDeg - oldDeg * oldDeg;

            if (newDeg > 0)
            {
                if (!p)
                    addKey(vertexKeys, vertexKeyIndex, v);
                vertexDegrees[v] = cast(size_t) newDeg;
            }
            else if (p)
            {
                vertexDegrees.remove(v);
                removeKey(vertexKeys, vertexKeyIndex, v);
            }
        }
        nVertices = vertexDegrees.length;

        // Update edge degrees (6 edges per tet)
        foreach (i; 0 .. 4)
            foreach (j; i + 1 .. 4)
            {
                auto ek = edgeKey(verts[i], verts[j]);
                auto p = ek in edgeDegrees;
                long oldDeg = p ? cast(long) *p : 0;
                long newDeg = oldDeg + delta;

                edgeDegreeSum += delta;
                edgeDegreeSqSum += newDeg * newDeg - oldDeg * oldDeg;

                if (newDeg > 0)
                {
                    if (!p)
                        addKey(edgeKeys, edgeKeyIndex, ek);
                    edgeDegrees[ek] = cast(size_t) newDeg;
                }
                else if (p)
                {
                    edgeDegrees.remove(ek);
                    removeKey(edgeKeys, edgeKeyIndex, ek);
                }
            }
        nEdges = edgeDegrees.length;
    }

    /// Variance of vertex degrees (number of tets per vertex)
    double vertexDegreeVariance() const
    {
        if (nVertices == 0) return 0.0;
        double n = cast(double) nVertices;
        double mean = cast(double) vertexDegreeSum / n;
        return cast(double) vertexDegreeSqSum / n - mean * mean;
    }

    /// Variance of edge degrees (number of tets per edge)
    double edgeDegreeVariance() const
    {
        if (nEdges == 0) return 0.0;
        double n = cast(double) nEdges;
        double mean = cast(double) edgeDegreeSum / n;
        return cast(double) edgeDegreeSqSum / n - mean * mean;
    }

    /// Recompute all counts and degree maps from scratch.
    void recount()
    {
        nTetrahedra = tets.length;

        // Clear all index structures
        vertexDegrees = typeof(vertexDegrees).init;
        edgeDegrees = typeof(edgeDegrees).init;
        vertexKeys = null;
        vertexKeyIndex = typeof(vertexKeyIndex).init;
        edgeKeys = null;
        edgeKeyIndex = typeof(edgeKeyIndex).init;
        vertexDegreeSum = 0;
        vertexDegreeSqSum = 0;
        edgeDegreeSum = 0;
        edgeDegreeSqSum = 0;

        // Clear face index
        faceTets = typeof(faceTets).init;
        faceKeys = null;
        faceKeyIndex = typeof(faceKeyIndex).init;

        // Rebuild face index
        foreach (i, ref t; tets)
            foreach (fi; 0 .. 4)
            {
                auto fk = t.face(cast(ubyte) fi);
                if (fk !in faceTets)
                    addKey(faceKeys, faceKeyIndex, fk);
                faceTets[fk] ~= i;
            }

        // Rebuild degree maps (also rebuilds vertexKeys/edgeKeys)
        foreach (ref t; tets)
            adjustDegrees(t.vertices, +1);

        // nVertices and nEdges are set by adjustDegrees
        nTriangles = faceTets.length;

        // Update next vertex counter
        static if (__traits(isIntegral, VertexLabel))
        {
            nextVertexCounter = 0;
            foreach (v; vertexDegrees.byKey())
                if (v >= nextVertexCounter)
                    nextVertexCounter = cast(VertexLabel)(v + 1);
        }
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
    /// the given vertex labels. Uses indexes for O(1) lookups.
    bool hasSimplexContaining(const VertexLabel[] vertices, const size_t[] excludeIndices = null) const
    {
        if (vertices.length == 1)
        {
            if (excludeIndices is null)
                return (vertices[0] in vertexDegrees) !is null;
            // Check if vertex exists in any non-excluded tet
            auto p = vertices[0] in vertexDegrees;
            if (p is null) return false;
            // Need to check if all tets containing this vertex are excluded
            // Fall through to linear scan for this rare case
        }
        else if (vertices.length == 2)
        {
            auto ek = edgeKey(vertices[0], vertices[1]);
            auto p = ek in edgeDegrees;
            if (p is null) return false;
            if (excludeIndices is null) return true;
            // Count how many excluded tets contain this edge
            size_t excludedCount = 0;
            foreach (idx; excludeIndices)
                if (idx < tets.length && canFind(tets[idx].vertices[], vertices[0])
                    && canFind(tets[idx].vertices[], vertices[1]))
                    excludedCount++;
            return *p > excludedCount;
        }
        else if (vertices.length == 3)
        {
            VertexLabel[3] fk = [vertices[0], vertices[1], vertices[2]];
            sort(fk[]);
            auto p = fk in faceTets;
            if (p is null) return false;
            if (excludeIndices is null) return (*p).length > 0;
            foreach (tetIdx; *p)
                if (!canFind(excludeIndices, tetIdx))
                    return true;
            return false;
        }

        // Fallback: linear scan
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

    /// Remove tetrahedron at the given index (swaps with last element).
    /// Updates all indexes (face, degree maps, key arrays).
    void removeTetrahedron(size_t idx)
    in (idx < tets.length)
    {
        size_t lastIdx = tets.length - 1;

        // Update degree maps for removed tet
        adjustDegrees(tets[idx].vertices, -1);

        // Step 1: Remove idx from face index for the tet being deleted
        foreach (fi; 0 .. 4)
        {
            auto fk = tets[idx].face(cast(ubyte) fi);
            auto p = fk in faceTets;
            if (p)
            {
                auto pos = countUntil(*p, idx);
                if (pos >= 0)
                {
                    *p = (*p).remove(cast(size_t) pos);
                    if ((*p).length == 0)
                    {
                        faceTets.remove(fk);
                        removeKey(faceKeys, faceKeyIndex, fk);
                    }
                }
            }
        }

        // Step 2: If not removing the last element, update lastIdx → idx in face index
        if (idx != lastIdx)
        {
            foreach (fi; 0 .. 4)
            {
                auto fk = tets[lastIdx].face(cast(ubyte) fi);
                auto p = fk in faceTets;
                if (p)
                {
                    foreach (ref entry; *p)
                        if (entry == lastIdx)
                            entry = idx;
                }
            }
        }

        // Step 3: Swap and shrink
        tets[idx] = tets[lastIdx];
        tets = tets[0 .. lastIdx];

        nTriangles = faceTets.length;
        nTetrahedra = tets.length;
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

        // Check that the replacement tet [outer] doesn't already exist as a surviving tet
        if (hasSimplexContaining(outer, tetIndices))
            return false;

        // Remove the 4 tets (remove from highest index first to preserve indices)
        sort!"a > b"(tetIndices);
        foreach (idx; tetIndices)
            removeTetrahedron(idx);

        // Add the single replacement tetrahedron
        VertexLabel[4] newVerts = outer[0 .. 4];
        addTetrahedron(newVerts);

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
        VertexLabel[4] newTet1 = [outer[0], outer[1], outer[2], d];
        VertexLabel[4] newTet2 = [outer[0], outer[1], outer[2], e];
        addTetrahedron(newTet1);
        addTetrahedron(newTet2);

        return true;
    }

    /// Returns a fresh vertex label not currently used. O(1).
    /// Only available for integral VertexLabel types.
    static if (__traits(isIntegral, VertexLabel))
    VertexLabel nextVertexLabel()
    {
        return nextVertexCounter++;
    }

    /// Collect all distinct vertex labels used in the triangulation
    VertexLabel[] vertexLabels() const
    {
        auto labels = vertexDegrees.keys;
        sort(labels);
        return labels;
    }

    /// Collect all distinct edges (pairs of vertex labels) in the triangulation
    VertexLabel[2][] edges() const
    {
        return edgeDegrees.keys;
    }

    /**
     * Returns all valid Pachner moves for the current triangulation.
     * Uses face/edge/vertex indexes for efficiency.
     */
    PachnerMove!VertexLabel[] validMoves() const
    {
        PachnerMove!VertexLabel[] moves;

        // 1-4 moves: one per tetrahedron (always valid since we use a fresh vertex)
        static if (__traits(isIntegral, VertexLabel))
        {
            VertexLabel fresh = nextVertexCounter;
            foreach (i; 0 .. tets.length)
                moves ~= PachnerMove!VertexLabel.make14(i, fresh);
        }

        // 2-3 moves: iterate faces with exactly 2 tets
        foreach (ref kv; faceTets.byKeyValue())
        {
            if (kv.value.length != 2)
                continue;
            auto i = kv.value[0];
            auto j = kv.value[1];
            auto face = kv.key;

            // Find non-shared vertex from each tet
            VertexLabel d, e;
            foreach (u; tets[i].vertices)
                if (!canFind(face[], u)) { d = u; break; }
            foreach (u; tets[j].vertices)
                if (!canFind(face[], u)) { e = u; break; }

            // Check validity: edge [d,e] must not exist in any tet
            // (tets i,j have vertices [a,b,c,d] and [a,b,c,e], neither contains both d and e)
            auto ek = edgeKey(d, e);
            if (ek !in edgeDegrees)
                moves ~= PachnerMove!VertexLabel.make23(i, j, face, d, e);
        }

        // 3-2 moves: iterate edges with degree exactly 3
        foreach (ref kv; edgeDegrees.byKeyValue())
        {
            if (kv.value != 3)
                continue;
            auto edge = kv.key;

            // Find the 3 tets containing this edge
            size_t[3] tetIndices;
            size_t count = 0;
            foreach (i, ref t; tets)
            {
                if (canFind(t.vertices[], edge[0]) && canFind(t.vertices[], edge[1]))
                {
                    if (count < 3) tetIndices[count] = i;
                    count++;
                }
            }
            if (count != 3) continue;

            // Collect outer vertices
            VertexLabel[3] outer;
            size_t outerCount = 0;
            foreach (idx; tetIndices)
                foreach (u; tets[idx].vertices)
                    if (u != edge[0] && u != edge[1] && !canFind(outer[0 .. outerCount], u))
                        outer[outerCount++] = u;

            if (outerCount != 3) continue;

            // Check validity: face [a,b,c] must not exist in surviving tets
            VertexLabel[3] fk = outer;
            sort(fk[]);
            auto fp = fk in faceTets;
            if (fp is null)
            {
                moves ~= PachnerMove!VertexLabel.make32(edge, outer);
            }
            else
            {
                // Check if all tets with this face are among the 3 being removed
                bool survivorHasFace = false;
                foreach (tetIdx; *fp)
                    if (tetIdx != tetIndices[0] && tetIdx != tetIndices[1] && tetIdx != tetIndices[2])
                        { survivorHasFace = true; break; }
                if (!survivorHasFace)
                    moves ~= PachnerMove!VertexLabel.make32(edge, outer);
            }
        }

        // 4-1 moves: iterate vertices with degree exactly 4
        foreach (ref kv; vertexDegrees.byKeyValue())
        {
            if (kv.value != 4)
                continue;
            auto v = kv.key;

            // Find the 4 tets containing this vertex
            size_t[4] tetIndices;
            size_t count = 0;
            foreach (i, ref t; tets)
                if (canFind(t.vertices[], v))
                {
                    if (count < 4) tetIndices[count] = i;
                    count++;
                }
            if (count != 4) continue;

            // Collect outer vertices
            VertexLabel[4] outer;
            size_t outerCount = 0;
            foreach (idx; tetIndices)
                foreach (u; tets[idx].vertices)
                    if (u != v && !canFind(outer[0 .. outerCount], u))
                        outer[outerCount++] = u;

            if (outerCount != 4) continue;

            // Verify each tet has exactly 3 outer vertices
            bool valid = true;
            foreach (idx; tetIndices)
            {
                int oc = 0;
                foreach (u; tets[idx].vertices)
                    if (canFind(outer[], u))
                        oc++;
                if (oc != 3) { valid = false; break; }
            }
            if (!valid) continue;

            // Check that replacement tet [outer] doesn't already exist as a surviving tet
            if (hasSimplexContaining(outer[], tetIndices[]))
                continue;

            moves ~= PachnerMove!VertexLabel.make41(v, outer);
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
        copy.vertexDegrees = cast(size_t[VertexLabel]) vertexDegrees.dup;
        copy.edgeDegrees = cast(size_t[VertexLabel[2]]) edgeDegrees.dup;
        copy.vertexDegreeSum = vertexDegreeSum;
        copy.vertexDegreeSqSum = vertexDegreeSqSum;
        copy.edgeDegreeSum = edgeDegreeSum;
        copy.edgeDegreeSqSum = edgeDegreeSqSum;

        // Deep copy face index
        foreach (ref kv; faceTets.byKeyValue())
            copy.faceTets[kv.key] = kv.value.dup;

        // Copy parallel key arrays
        copy.vertexKeys = vertexKeys.dup;
        copy.vertexKeyIndex = cast(size_t[VertexLabel]) vertexKeyIndex.dup;
        copy.edgeKeys = edgeKeys.dup;
        copy.edgeKeyIndex = cast(size_t[VertexLabel[2]]) edgeKeyIndex.dup;
        copy.faceKeys = faceKeys.dup;
        copy.faceKeyIndex = cast(size_t[VertexLabel[3]]) faceKeyIndex.dup;

        static if (__traits(isIntegral, VertexLabel))
            copy.nextVertexCounter = nextVertexCounter;

        return copy;
    }

    /**
     * Propose a random Pachner move in O(1) expected time.
     *
     * Picks a move type weighted by candidate pool sizes, then picks
     * a random candidate. Returns null if the candidate is invalid
     * (caller should treat as a rejected step in Metropolis).
     *
     * Pool sizes: 1-4 = nTets, 2-3 = nFaces, 3-2 = nEdges, 4-1 = nVertices.
     * Within each pool, a random element is selected and checked for validity.
     */
    Nullable!(PachnerMove!VertexLabel) proposeMove(Rng)(ref Rng rng)
    {
        import std.random : uniform;
        alias PM = PachnerMove!VertexLabel;

        static if (__traits(isIntegral, VertexLabel))
            size_t n14 = tets.length;
        else
            size_t n14 = 0;

        size_t n23 = faceKeys.length;
        size_t n32 = edgeKeys.length;
        size_t n41 = vertexKeys.length;
        size_t total = n14 + n23 + n32 + n41;

        if (total == 0)
            return typeof(return).init;

        auto pick = uniform(0, total, rng);

        if (pick < n14)
        {
            // 1-4 move: pick random tet (always valid)
            static if (__traits(isIntegral, VertexLabel))
            {
                auto tetIdx = uniform(0, tets.length, rng);
                return nullable(PM.make14(tetIdx, nextVertexLabel()));
            }
        }
        else if (pick < n14 + n23)
        {
            // 2-3 move: pick random face, check if it has exactly 2 tets
            auto fk = faceKeys[uniform(0, faceKeys.length, rng)];
            auto p = fk in faceTets;
            if (p is null || (*p).length != 2)
                return typeof(return).init;

            auto i = (*p)[0];
            auto j = (*p)[1];

            // Find non-shared vertices
            VertexLabel d, e;
            foreach (u; tets[i].vertices)
                if (!canFind(fk[], u)) { d = u; break; }
            foreach (u; tets[j].vertices)
                if (!canFind(fk[], u)) { e = u; break; }

            // Check validity: edge [d,e] must not exist
            auto ek = edgeKey(d, e);
            if (ek in edgeDegrees)
                return typeof(return).init;

            return nullable(PM.make23(i, j, fk, d, e));
        }
        else if (pick < n14 + n23 + n32)
        {
            // 3-2 move: pick random edge, check degree == 3
            auto edge = edgeKeys[uniform(0, edgeKeys.length, rng)];
            auto degP = edge in edgeDegrees;
            if (degP is null || *degP != 3)
                return typeof(return).init;

            // Find the 3 tets containing this edge
            size_t[3] tetIndices;
            size_t count = 0;
            foreach (i, ref t; tets)
                if (canFind(t.vertices[], edge[0]) && canFind(t.vertices[], edge[1]))
                {
                    if (count < 3) tetIndices[count] = i;
                    count++;
                }
            if (count != 3) return typeof(return).init;

            // Collect outer vertices
            VertexLabel[3] outer;
            size_t outerCount = 0;
            foreach (idx; tetIndices)
                foreach (u; tets[idx].vertices)
                    if (u != edge[0] && u != edge[1] && !canFind(outer[0 .. outerCount], u))
                        outer[outerCount++] = u;
            if (outerCount != 3) return typeof(return).init;

            // Check validity: face [a,b,c] must not exist in surviving tets
            VertexLabel[3] fk = outer;
            sort(fk[]);
            auto fp = fk in faceTets;
            if (fp !is null)
            {
                foreach (tetIdx; *fp)
                    if (tetIdx != tetIndices[0] && tetIdx != tetIndices[1] && tetIdx != tetIndices[2])
                        return typeof(return).init;
            }

            return nullable(PM.make32(edge, outer));
        }
        else
        {
            // 4-1 move: pick random vertex, check degree == 4
            auto v = vertexKeys[uniform(0, vertexKeys.length, rng)];
            auto degP = v in vertexDegrees;
            if (degP is null || *degP != 4)
                return typeof(return).init;

            // Find the 4 tets containing this vertex
            size_t[4] tetIndices;
            size_t count = 0;
            foreach (i, ref t; tets)
                if (canFind(t.vertices[], v))
                {
                    if (count < 4) tetIndices[count] = i;
                    count++;
                }
            if (count != 4) return typeof(return).init;

            // Collect outer vertices
            VertexLabel[4] outer;
            size_t outerCount = 0;
            foreach (idx; tetIndices)
                foreach (u; tets[idx].vertices)
                    if (u != v && !canFind(outer[0 .. outerCount], u))
                        outer[outerCount++] = u;
            if (outerCount != 4) return typeof(return).init;

            // Verify each tet has exactly 3 outer vertices
            foreach (idx; tetIndices)
            {
                int oc = 0;
                foreach (u; tets[idx].vertices)
                    if (canFind(outer[], u))
                        oc++;
                if (oc != 3) return typeof(return).init;
            }

            // Check that replacement tet [outer] doesn't already exist
            if (hasSimplexContaining(outer[], tetIndices[]))
                return typeof(return).init;

            return nullable(PM.make41(v, outer));
        }

        return typeof(return).init;
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
    // Degree tracking for 4-simplex boundary
    auto tri = fourSimplexBoundary();

    // Each vertex is in 4 tets (every tet except the one that omits it)
    foreach (v; 0 .. 5)
        assert(tri.vertexDegrees[v] == 4);

    // Each edge is in 3 tets (every tet containing both endpoints)
    assert(tri.edgeDegrees.length == 10);
    foreach (ref pair; tri.edgeDegrees.byKeyValue())
        assert(pair.value == 3);

    // Variance should be 0 since all degrees are equal
    import std.math : abs;
    assert(abs(tri.vertexDegreeVariance()) < 1e-10);
    assert(abs(tri.edgeDegreeVariance()) < 1e-10);
}

unittest
{
    // After 1-4 move, degrees become non-uniform
    auto tri = fourSimplexBoundary();
    tri.move14(0, 5);

    // New vertex 5 should be in exactly 4 tets
    assert(tri.vertexDegrees[5] == 4);

    // Variance should be > 0 since degrees are now non-uniform
    assert(tri.vertexDegreeVariance() > 0);
    assert(tri.edgeDegreeVariance() > 0);
}

unittest
{
    import std.random : Random, uniform;
    import std.math : abs;

    // Verify degree variance tracking matches recount through random walk
    auto tri = fourSimplexBoundary();
    auto rng = Random(99);
    alias Move = PachnerMove!size_t;

    foreach (step; 0 .. 50)
    {
        auto moves = tri.validMoves();
        auto idx = uniform(0, moves.length, rng);
        auto move = moves[idx];
        if (move.type == Move.Type.move14)
            move.newVertex = tri.nextVertexLabel();
        move.execute(tri);

        // Save tracked values
        auto vVar = tri.vertexDegreeVariance();
        auto eVar = tri.edgeDegreeVariance();
        auto vSum = tri.vertexDegreeSum;
        auto vSqSum = tri.vertexDegreeSqSum;
        auto eSum = tri.edgeDegreeSum;
        auto eSqSum = tri.edgeDegreeSqSum;

        // Recount and compare
        tri.recount();
        assert(tri.vertexDegreeSum == vSum);
        assert(tri.vertexDegreeSqSum == vSqSum);
        assert(tri.edgeDegreeSum == eSum);
        assert(tri.edgeDegreeSqSum == eSqSum);
        assert(abs(tri.vertexDegreeVariance() - vVar) < 1e-10);
        assert(abs(tri.edgeDegreeVariance() - eVar) < 1e-10);
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
        assert(ok, "Move execution failed");
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

unittest
{
    import std.random : Random, uniform;

    // Test proposeMove: every non-null proposal must be executable
    // and must correspond to a valid move.
    auto tri = fourSimplexBoundary();
    auto rng = Random(77);
    alias Move = PachnerMove!size_t;

    size_t proposals = 0;
    size_t accepted = 0;

    foreach (step; 0 .. 500)
    {
        // Get all valid moves for comparison
        auto allValid = tri.validMoves();
        assert(allValid.length > 0);

        // Try proposing
        auto proposed = tri.proposeMove(rng);
        proposals++;

        if (!proposed.isNull)
        {
            auto move = proposed.get;

            // Verify the move executes successfully on a copy
            auto copy = tri.dup;
            if (move.type == Move.Type.move14)
                move.newVertex = copy.nextVertexLabel();
            bool ok = move.execute(copy);
            assert(ok, "Proposed move failed to execute");

            // Execute on the real triangulation
            if (move.type == Move.Type.move14)
                move.newVertex = tri.nextVertexLabel();
            ok = move.execute(tri);
            assert(ok);
            accepted++;

            // Verify consistency
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
        else
        {
            // Null proposal = rejected step, pick from validMoves instead
            auto idx = uniform(0, allValid.length, rng);
            auto move = allValid[idx];
            if (move.type == Move.Type.move14)
                move.newVertex = tri.nextVertexLabel();
            move.execute(tri);
        }
    }

    // Should have accepted a reasonable fraction of proposals
    assert(accepted > 50, "Too few proposals accepted");
}

unittest
{
    import std.random : Random, uniform;
    import std.math : abs;

    // Test face index consistency through random walk
    auto tri = fourSimplexBoundary();
    auto rng = Random(55);
    alias Move = PachnerMove!size_t;

    foreach (step; 0 .. 200)
    {
        auto moves = tri.validMoves();
        auto idx = uniform(0, moves.length, rng);
        auto move = moves[idx];
        if (move.type == Move.Type.move14)
            move.newVertex = tri.nextVertexLabel();
        move.execute(tri);

        // Verify face index matches brute-force
        size_t[][size_t[3]] expectedFaces;
        foreach (i, ref t; tri.tets)
            foreach (fi; 0 .. 4)
            {
                auto fk = t.face(cast(ubyte) fi);
                expectedFaces[fk] ~= i;
            }

        assert(tri.faceTets.length == expectedFaces.length,
            "Face index size mismatch");
        assert(tri.nTriangles == expectedFaces.length,
            "nTriangles mismatch");
        assert(tri.faceKeys.length == expectedFaces.length,
            "faceKeys size mismatch");

        // Verify degree tracking
        auto vVar = tri.vertexDegreeVariance();
        auto eVar = tri.edgeDegreeVariance();
        tri.recount();
        assert(abs(tri.vertexDegreeVariance() - vVar) < 1e-10);
        assert(abs(tri.edgeDegreeVariance() - eVar) < 1e-10);
    }
}
