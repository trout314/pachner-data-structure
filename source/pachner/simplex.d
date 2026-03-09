/**
 * Simplex representation for 3-manifold triangulations.
 *
 * A 3-manifold triangulation is built from tetrahedra (3-simplices).
 * Each tetrahedron has 4 vertices. Gluings between tetrahedra are
 * determined implicitly by shared vertex labels — two tetrahedra
 * sharing 3 vertex labels are glued along that common face.
 */
module pachner.simplex;

import std.algorithm : sort;

/**
 * A tetrahedron (3-simplex) defined by 4 vertex labels.
 *
 * Templated on the vertex label type for flexibility (e.g., size_t,
 * int, string, or a custom identifier type).
 */
struct Tetrahedron(VertexLabel = size_t)
{
    /// The 4 vertex labels of this tetrahedron
    VertexLabel[4] vertices;

    /// Returns the face (as sorted vertex labels) opposite to local vertex v
    VertexLabel[3] face(ubyte v) const
    in (v < 4)
    {
        VertexLabel[3] f;
        size_t k = 0;
        foreach (i; 0 .. 4)
            if (i != v)
                f[k++] = vertices[i];
        sort(f[]);
        return f;
    }

    /// Returns a copy of the vertex labels in sorted order
    VertexLabel[4] sortedVertices() const
    {
        VertexLabel[4] s = vertices;
        sort(s[]);
        return s;
    }

    /// Compare by sorted vertex labels (order-independent)
    int opCmp(ref const Tetrahedron rhs) const
    {
        auto a = sortedVertices();
        auto b = rhs.sortedVertices();
        foreach (i; 0 .. 4)
        {
            if (a[i] < b[i]) return -1;
            if (a[i] > b[i]) return 1;
        }
        return 0;
    }

    /// Returns all 4 faces (each as sorted vertex labels)
    VertexLabel[3][4] faces() const
    {
        VertexLabel[3][4] result;
        foreach (v; 0 .. 4)
            result[v] = face(cast(ubyte) v);
        return result;
    }
}

unittest
{
    auto t = Tetrahedron!size_t([0, 1, 2, 3]);

    auto f = t.face(0);
    assert(f == [1, 2, 3]);

    f = t.face(3);
    assert(f == [0, 1, 2]);

    auto allFaces = t.faces();
    assert(allFaces[0] == [1, 2, 3]);
    assert(allFaces[1] == [0, 2, 3]);
    assert(allFaces[2] == [0, 1, 3]);
    assert(allFaces[3] == [0, 1, 2]);
}

unittest
{
    // Test with string vertex labels
    auto t = Tetrahedron!string(["a", "b", "c", "d"]);

    auto f = t.face(0);
    assert(f == ["b", "c", "d"]);
}
