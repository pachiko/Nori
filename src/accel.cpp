/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <numeric>
#include <array>
#include <chrono>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    if (m_mesh == nullptr) return;
    uint32_t count = m_mesh->getTriangleCount();
    std::vector<uint32_t> tris(count);
    std::iota(tris.begin(), tris.end(), 0);
    depth = numInterior = numLeaf = numTris = 0;

    auto start = std::chrono::high_resolution_clock::now();
    m_root = recursiveBuild(m_bbox, tris);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "OctTree build time: " << duration.count() << " (ms) \n";
    std::cout << "OctTree number of interior nodes: " << numInterior << "\n";
    std::cout << "OctTree number of leaf nodes: " << numLeaf << "\n";
    std::cout << "OctTree average number of triangles per leaf node: " << numTris/numLeaf << "\n";
}

std::unique_ptr<OctTreeNode> Accel::recursiveBuild(const BoundingBox3f& bbox, std::vector<uint32_t>& tris) {
    if (tris.empty()) return nullptr;
    depth++;
    auto res = std::make_unique<OctTreeNode>(bbox);

    if (tris.size() <= NORI_NODE_MAX_TRI_COUNT || depth >= NORI_NODE_MAX_TREE_DEPTH) {
        res->triIndices = std::move(tris);
        depth--;
        numLeaf++;
        numTris += res->triIndices.size();
        return res;
    }

    std::vector<uint32_t> tri_list[8];
    std::vector<BoundingBox3f> childBoxes;
    Point3f center = bbox.getCenter();

    for (int i = 0; i < 8; i++) {
        Point3f corner = bbox.getCorner(i);

        Point3f min(std::min(corner[0], center[0]),
                    std::min(corner[1], center[1]),
                    std::min(corner[2], center[2]));

        Point3f max(std::max(corner[0], center[0]),
            std::max(corner[1], center[1]),
            std::max(corner[2], center[2]));

        childBoxes.push_back(BoundingBox3f(min, max));
    }

    for (auto tri : tris) {
        const BoundingBox3f& triB = m_mesh->getBoundingBox(tri);
        for (int i = 0; i < 8; ++i) {
            const BoundingBox3f& childB = childBoxes[i];
            if (triB.overlaps(childB)) tri_list[i].push_back(tri);
        }
    } 

    res->children.resize(8);
    for (int i = 0; i < 8; ++i) res->children[i] = recursiveBuild(childBoxes[i], tri_list[i]);
    depth--;
    numInterior++;
    return res;
}


void OctTreeNode::rayIntersect(const Mesh& mesh, Ray3f& ray, Intersection& its, bool& hit, uint32_t& triIdx, bool shadowRay) {
    bool hitBound;
    float nearT;
    float farT;

    hitBound = bound.rayIntersect(ray, nearT, farT);
    if (!hitBound || nearT > ray.maxt) return;

    if (!children.empty()) {
        std::array<float, 8> times;
        for (int i = 0; i < 8; i++) {
            if (children[i]) {
                hitBound = children[i]->bound.rayIntersect(ray, nearT, farT);
                if (hitBound && ray.maxt > nearT) {
                    times[i] = nearT;
                    continue;
                }
            }
            times[i] = std::numeric_limits<float>::infinity();
        }

        std::array<int, 8> idx;
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&times](int i1, int i2) { return times[i1] < times[i2]; });

        for (auto i : idx) {
            if (children[i] && (!hit || times[i] < ray.maxt)) {
                children[i]->rayIntersect(mesh, ray, its, hit, triIdx, shadowRay);
                if (shadowRay && hit) return;
            }
        }
    }
    else {
        for (auto tri : triIndices) {
            float u, v, t;
            if (mesh.rayIntersect(tri, ray, u, v, t)) {
                hit = true;
                if (shadowRay) return;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = &mesh;
                triIdx = tri;
            }
        }
    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    m_root->rayIntersect(*m_mesh, ray, its, foundIntersection, f, shadowRay);

    if (!shadowRay && foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;
        its.bary = bary;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
        its.tri_index = Point3f(idx0, idx1, idx2);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

