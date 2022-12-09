/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once
#include <memory>
#include <nori/mesh.h>
#include <tbb/tbb.h>
#include <atomic>

#define NORI_NODE_MAX_TRI_COUNT 10 /* Maximum number of triangles in a node */
#define NORI_NODE_MAX_TREE_DEPTH 10 /* Maximum depth of the oct tree*/

NORI_NAMESPACE_BEGIN

/**
 * \brief Node data structure
 *
 * Contains 8 children nodes in a pointer OR the indices of triangles.
 * Cannot contain children nodes AND triangles.
 * Also contains the bounding box of the node
 */
struct OctTreeNode {
    std::vector<std::unique_ptr<OctTreeNode>> children;
    std::vector<uint32_t> triIndices;
    BoundingBox3f bound;

    OctTreeNode(const BoundingBox3f& b) : bound(b) { }
    void rayIntersect(const Mesh& mesh, Ray3f& ray, Intersection& its, bool& hit, uint32_t& triIdx, bool shadowRay);
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure
    void build();
    
    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    std::unique_ptr<OctTreeNode> m_root; // root node of the oct tree
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    
    std::atomic<uint32_t> numInterior;
    std::atomic<uint32_t> numLeaf;
    std::atomic<uint32_t> numTris;
    bool parallelBuildMode = true;

    std::unique_ptr<OctTreeNode> recursiveBuild(const BoundingBox3f& bbox, std::vector<uint32_t>& tris, uint8_t depth);
};

NORI_NAMESPACE_END
