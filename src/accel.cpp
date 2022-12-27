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
    m_meshes.push_back(mesh);
    m_bbox.expandBy(mesh->getBoundingBox());
}

void Accel::build() {
    if (m_meshes.empty()) return;

    uint32_t count = 0;
    for (auto mesh : m_meshes) count += mesh->getTriangleCount();

    numInterior = numLeaf = numTris = 0;

    std::vector<triMeshPair> tris;
    tris.reserve(count);
    for (uint32_t m = 0; m < m_meshes.size(); m++) {
        Mesh* mesh = m_meshes[m];
        uint32_t c = mesh->getTriangleCount();
        for (uint32_t i = 0; i < c; i++) {
            tris.push_back(triMeshPair(m, i));
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    m_root = recursiveBuild(m_bbox, tris, 0);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "OctTree build time: " << duration.count() << " (ms) \n";
    std::cout << "OctTree number of interior nodes: " << numInterior << "\n";
    std::cout << "OctTree number of leaf nodes: " << numLeaf << "\n";
    std::cout << "OctTree average number of triangles per leaf node: " << static_cast<double>(numTris)/numLeaf << "\n";
}

std::unique_ptr<OctTreeNode> Accel::recursiveBuild(const BoundingBox3f& bbox, std::vector<triMeshPair>& tris, uint8_t depth) {
    if (tris.empty()) return nullptr;
    auto res = std::make_unique<OctTreeNode>(bbox);

    if (tris.size() <= NORI_NODE_MAX_TRI_COUNT || depth >= NORI_NODE_MAX_TREE_DEPTH) {
        res->triIndices = std::move(tris);
        numLeaf++;
        numTris += res->triIndices.size();
        return res;
    }

    std::vector<triMeshPair> tri_list[8];
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
        const BoundingBox3f& triB = m_meshes[std::get<0>(tri)]->getBoundingBox(std::get<1>(tri));
        for (int i = 0; i < 8; ++i) {
            const BoundingBox3f& childB = childBoxes[i];
            if (triB.overlaps(childB)) tri_list[i].push_back(tri);
        }
    } 
    tris.clear();

    res->children.resize(8);

    if (parallelBuildMode) {
        tbb::task_group tg;
        for (int i = 0; i < 8; ++i) {
            tg.run([&, i]() {
                res->children[i] = recursiveBuild(childBoxes[i], tri_list[i], depth + 1);
            });
        }
        tg.wait();
    }
    else {
        for (int i = 0; i < 8; ++i) res->children[i] = recursiveBuild(childBoxes[i], tri_list[i], depth + 1);
    }

    numInterior++;
    return res;
}

void OctTreeNode::rayIntersect(const std::vector<Mesh*>& meshes, Ray3f& ray, Intersection& its, bool& hit,
    uint32_t& triIdx, bool shadowRay) {
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
                children[i]->rayIntersect(meshes, ray, its, hit, triIdx, shadowRay);
                if (shadowRay && hit) return;
            }
        }
    }
    else {
        for (auto tri : triIndices) {
            float u, v, t;
            Mesh* mesh = meshes[std::get<0>(tri)];
            if (mesh->rayIntersect(std::get<1>(tri), ray, u, v, t)) {
                hit = true;
                if (shadowRay) return;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = mesh;
                triIdx = std::get<1>(tri);
            }
        }
    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    m_root->rayIntersect(m_meshes, ray, its, foundIntersection, f, shadowRay);

    if (!shadowRay && foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

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

