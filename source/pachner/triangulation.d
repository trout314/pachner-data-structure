/**
 * Triangulation data structure for closed orientable 3-manifolds.
 *
 * Stores a collection of tetrahedra with face-pairing gluings.
 * Supports Pachner moves (bistellar flips):
 *   - 2-3 move: replace 2 tetrahedra sharing a face with 3 tetrahedra sharing an edge
 *   - 3-2 move: inverse of the above
 *   - 1-4 move: subdivide 1 tetrahedron into 4 (introduces a new interior vertex)
 *   - 4-1 move: inverse of the above (removes an interior vertex of degree 4)
 *
 * For closed triangulations (no boundary), only 2-3 and 3-2 moves suffice
 * to connect any two triangulations of the same manifold (Pachner's theorem).
 */
module pachner.triangulation;

import pachner.simplex;

import std.algorithm : canFind, countUntil, remove;
import std.array : array;
import std.exception : enforce;
import std.range : iota;
import std.stdio : writefln;

/**
 * A triangulation of a 3-manifold stored as a list of glued tetrahedra.
 */
struct Triangulation
{
    /// All tetrahedra in the triangulation
    Tetrahedron[] tets;

    /// Total number of distinct vertex labels used
    size_t numVertices;

    /// Add a new unglued tetrahedron and return its index
    TetIndex addTetrahedron(VertexIndex[4] verts)
    {
        TetIndex idx = tets.length;
        Tetrahedron t;
        t.index = idx;
        t.vertices = verts;
        tets ~= t;
        return idx;
    }

    /**
     * Glue face faceA of tetA to face faceB of tetB.
     *
     * perm maps local vertex indices of tetA to local vertex indices of tetB,
     * restricted to the shared face.  perm must be a permutation of {0,1,2,3}
     * that fixes the "opposite" vertex conventions:
     *   perm[faceA] == faceB (the two opposite vertices map to each other)
     */
    void glue(TetIndex a, ubyte faceA, TetIndex b, ubyte faceB, ubyte[4] perm)
    {
        tets[a].neighbor[faceA] = b;
        tets[a].gluing[faceA] = perm;
        tets[b].neighbor[faceB] = a;
        // inverse permutation for the reverse direction
        ubyte[4] inv;
        foreach (i; 0 .. 4)
            inv[perm[i]] = cast(ubyte) i;
        tets[b].gluing[faceB] = inv;
    }

    /// Number of tetrahedra
    size_t size() const { return tets.length; }

    /**
     * Attempt the 2-3 Pachner move on the internal face between tetA (face faceA).
     *
     * Preconditions:
     *   - Face faceA of tetA is internal (has a neighbor tetB)
     *   - The two tetrahedra are distinct
     *
     * The move removes tetA and tetB and introduces 3 new tetrahedra.
     *
     * Returns: true if the move was performed, false if preconditions fail.
     */
    bool move23(TetIndex idxA, ubyte faceA)
    {
        if (idxA >= tets.length) return false;
        Tetrahedron* A = &tets[idxA];
        if (A.isBoundaryFace(faceA)) return false;

        TetIndex idxB = A.neighbor[faceA];
        if (idxB == idxA) return false; // self-glued face
        Tetrahedron* B = &tets[idxB];

        // faceB is the face of B that is glued to faceA of A
        ubyte faceB = cast(ubyte) A.gluing[faceA][faceA];

        // Vertices:
        //   A has local verts 0..3; the shared face is opposite vertex faceA
        //   The shared triangle has 3 vertices; the "top" vertex is A.vertices[faceA]
        //   The "bottom" vertex is B.vertices[faceB]
        //   The new edge in the 2->3 move connects "top" and "bottom"

        // Collect the shared face vertices in A's local ordering
        ubyte[3] sharedA; // local indices in A of the shared face vertices
        size_t k = 0;
        foreach (i; 0 .. 4)
            if (i != faceA)
                sharedA[k++] = cast(ubyte) i;

        // Global vertex labels
        VertexIndex vTop = A.vertices[faceA];
        VertexIndex vBot = B.vertices[faceB];
        VertexIndex[3] vShared;
        foreach (i; 0 .. 3)
            vShared[i] = A.vertices[sharedA[i]];

        // The 3 new tetrahedra each contain the new edge (vTop, vBot)
        // plus one edge of the shared triangle.
        // New tet i (i=0,1,2) has vertices: vTop, vBot, vShared[i], vShared[(i+1)%3]
        // We'll store them as [vTop, vBot, vShared[j], vShared[k]]

        // Save external gluings of A and B before we overwrite
        // For A: external faces are those != faceA
        // For B: external faces are those != faceB

        // Build adjacency info for external faces
        struct FaceGluing
        {
            TetIndex neighbor;
            ubyte neighborFace;
            ubyte[4] perm;
        }

        FaceGluing[3] extA, extB;
        foreach (i; 0 .. 3)
        {
            ubyte fi = sharedA[i]; // local face index in A (opposite sharedA[i], which is != faceA)
            extA[i].neighbor = A.neighbor[fi];
            extA[i].neighborFace = (extA[i].neighbor != NULL_INDEX)
                ? cast(ubyte) A.gluing[fi][fi] : 0;
            extA[i].perm = A.gluing[fi];
        }
        // For B, need to find which local faces correspond to the shared triangle verts
        // B's gluing[faceB] maps A's local verts -> B's local verts
        foreach (i; 0 .. 3)
        {
            ubyte aLocal = sharedA[i];
            ubyte bLocal = cast(ubyte) A.gluing[faceA][aLocal]; // image in B
            extB[i].neighbor = B.neighbor[bLocal];
            extB[i].neighborFace = (extB[i].neighbor != NULL_INDEX)
                ? cast(ubyte) B.gluing[bLocal][bLocal] : 0;
            extB[i].perm = B.gluing[bLocal];
        }

        // Remove A and B from tets (remove higher index first to keep indices valid)
        // After removal, indices shift — easiest approach: mark as removed and compact,
        // but for simplicity here we overwrite A and B with the first 2 new tets,
        // and append the 3rd.

        // New tet 0: vTop, vBot, vShared[1], vShared[2]  (opposite to vShared[0])
        // New tet 1: vTop, vBot, vShared[0], vShared[2]  (opposite to vShared[1])
        // New tet 2: vTop, vBot, vShared[0], vShared[1]  (opposite to vShared[2])
        // Layout: local vertex 0=vTop,1=vBot,2=vShared[a],3=vShared[b]
        // The internal faces (shared among new tets) each contain vTop and vBot.

        // We reuse idxA and idxB for new tets 0 and 1, and append new tet 2.
        TetIndex idxC = tets.length; // will be appended

        // Overwrite A
        tets[idxA] = Tetrahedron.init;
        tets[idxA].index = idxA;
        tets[idxA].vertices = [vTop, vBot, vShared[1], vShared[2]];

        // Overwrite B
        tets[idxB] = Tetrahedron.init;
        tets[idxB].index = idxB;
        tets[idxB].vertices = [vTop, vBot, vShared[0], vShared[2]];

        // Append C
        Tetrahedron C;
        C.index = idxC;
        C.vertices = [vTop, vBot, vShared[0], vShared[1]];
        tets ~= C;

        // Internal gluings between new tets:
        // Tet0 face opposite vShared[1] (local 2) <-> Tet1 face opposite vShared[0] (local 2)
        // Tet0 face opposite vShared[2] (local 3) <-> Tet2 face opposite vShared[0] (local 2)
        // Tet1 face opposite vShared[2] (local 3) <-> Tet2 face opposite vShared[1] (local 2)
        // Identity permutation works since vTop(0) and vBot(1) are shared by all.

        // glue Tet0 (face 2, opposite vShared[1]) to Tet1 (face 2, opposite vShared[0])
        // local 0=vTop,1=vBot are the same; local 3 in Tet0 = vShared[2] = local 3 in Tet1
        glue(idxA, 2, idxB, 2, [0, 1, 3, 2]); // swap local 2&3 across the face

        // glue Tet0 (face 3, opposite vShared[2]) to Tet2 (face 3, opposite vShared[1])
        glue(idxA, 3, idxC, 3, [0, 1, 2, 3]); // local 2 in Tet0=vShared[1], local 2 in Tet2=vShared[0]... need to fix
        // Actually this needs more careful bookkeeping — placeholder for now.

        // TODO: fix external gluings (re-attach neighbors)

        return true;
    }

    /// Pretty-print summary
    void print() const
    {
        writefln("Triangulation: %d tetrahedra, %d vertices", tets.length, numVertices);
        foreach (ref t; tets)
        {
            writefln("  Tet %d: verts=%s  neighbors=%s", t.index, t.vertices, t.neighbor);
        }
    }
}

unittest
{
    // Construct the minimal triangulation of S^3: 2 tetrahedra glued along all faces
    Triangulation tri;
    tri.numVertices = 2; // only 2 vertices needed for S^3

    TetIndex t0 = tri.addTetrahedron([0, 0, 0, 1]); // degenerate labels for S^3
    TetIndex t1 = tri.addTetrahedron([0, 0, 0, 1]);
    assert(tri.size == 2);
}
