/**
 * Simplex representation for 3-manifold triangulations.
 *
 * A 3-manifold triangulation is built from tetrahedra (3-simplices).
 * Each tetrahedron has 4 vertices, 6 edges, 4 triangular faces (2-simplices),
 * and is glued to neighbors across its faces.
 */
module pachner.simplex;

import std.algorithm : sort, canFind;
import std.array : array;
import std.exception : enforce;

/// A vertex label (non-negative integer index)
alias VertexIndex = size_t;

/// A tetrahedron index within the triangulation
alias TetIndex = size_t;

/// Invalid/null index sentinel
enum size_t NULL_INDEX = size_t.max;

/**
 * Represents a single tetrahedron (3-simplex) in a triangulation.
 *
 * Vertices are labeled 0..3 locally. Face i is opposite vertex i,
 * i.e., the face spanned by vertices {0,1,2,3} \ {i}.
 *
 * Gluings: neighbor[i] is the tetrahedron glued along face i,
 * and gluing[i] maps local vertex indices to neighbor[i]'s local indices.
 */
struct Tetrahedron
{
    /// Global index of this tetrahedron in the triangulation
    TetIndex index;

    /// Global vertex labels for the 4 local vertices (0..3)
    VertexIndex[4] vertices;

    /// Neighboring tetrahedron index across each face (NULL_INDEX if boundary)
    TetIndex[4] neighbor = [NULL_INDEX, NULL_INDEX, NULL_INDEX, NULL_INDEX];

    /**
     * Gluing permutation for each face.
     * gluing[i][j] = local vertex in neighbor[i] corresponding to local vertex j
     * of this tetrahedron (restricted to the shared face).
     * Stored as a permutation of {0,1,2,3}.
     */
    ubyte[4][4] gluing;

    /// Returns the face (as sorted vertex indices) opposite to local vertex v
    VertexIndex[3] face(ubyte v) const
    in (v < 4)
    {
        VertexIndex[3] f;
        size_t k = 0;
        foreach (i; 0 .. 4)
            if (i != v)
                f[k++] = vertices[i];
        sort(f[]);
        return f;
    }

    /// True if face i is a boundary face (no neighbor)
    bool isBoundaryFace(ubyte i) const
    in (i < 4)
    {
        return neighbor[i] == NULL_INDEX;
    }
}

unittest
{
    Tetrahedron t;
    t.index = 0;
    t.vertices = [0, 1, 2, 3];

    auto f = t.face(0);
    assert(f == [1, 2, 3]);

    f = t.face(3);
    assert(f == [0, 1, 2]);

    assert(t.isBoundaryFace(0));
    t.neighbor[0] = 1;
    assert(!t.isBoundaryFace(0));
}
