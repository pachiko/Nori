#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>

NORI_NAMESPACE_BEGIN


// Material Sampling Path Tracer
class PathMats : public Integrator {
public:
    PathMats(const PropertyList& props) { }

    /// Compute the radiance value for a given ray. Just return green here
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f res(0.f); // accumulated radiance along the path

        float eta = 1.f; // product of all eta terms along the path
        Color3f throughput(1.f); // Throughput = BRDF * cos / pRR / p
        Ray3f r = ray; // ray

        for (uint32_t bounces = 0; ; bounces++) {
            Intersection its;
            if (!scene->rayIntersect(r, its)) break;

            if (its.mesh->isEmitter()) {
                EmitterQueryRecord rec;
                rec.origin = r.o;
                rec.hit = its.p;
                rec.hitN = its.shFrame.n;
                rec.light = its.mesh;
                res += its.mesh->getEmitter()->eval(rec) * throughput;
            }

            BSDFQueryRecord rec(its.shFrame.toLocal(-r.d));
            const BSDF* f = its.mesh->getBSDF();
            Color3f F = f->sample(rec, sampler->next2D());
            eta *= rec.eta;
            throughput *= F;
            r = Ray3f(its.p, its.shFrame.toWorld(rec.wo));

            float pRR = 1.f;
            if (bounces > 3) pRR = std::min(throughput.maxCoeff() * eta * eta, 0.99f);
            if (sampler->next1D() >= pRR) break;
            throughput /= pRR;
        }

        return res;
    }

    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return tfm::format(
            "PathMats"
        );
    }
};

NORI_REGISTER_CLASS(PathMats, "path_mats");
NORI_NAMESPACE_END