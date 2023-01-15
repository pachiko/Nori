﻿#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>

NORI_NAMESPACE_BEGIN


// Emitter Sampling Path Tracer. Basically Next Event Estimation
class PathEms : public Integrator {
public:
    PathEms(const PropertyList& props) { }

    /// Compute the radiance value for a given ray. Just return green here
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f res(0.f); // accumulated radiance along the path

        float eta = 1.f; // product of all eta terms along the path
        Color3f throughput(1.f); // Throughput = BRDF * cos / pRR / p
        Ray3f r = ray; // ray
        bool specular = false; // Was the previous hit's BSDF a Dirac-delta function?

        for (uint32_t bounces = 0; ; bounces++) {
            Intersection its;
            if (!scene->rayIntersect(r, its)) break;


            const BSDF* f = its.mesh->getBSDF();

            if (its.mesh->isEmitter() && (specular || bounces == 0)) {
                EmitterQueryRecord rec;
                rec.origin = r.o;
                rec.hit = its.p;
                rec.hitN = its.shFrame.n;
                rec.light = its.mesh;
                Color3f Le = its.mesh->getEmitter()->eval(rec);
                if (Le.isValid()) res += Le * throughput;
            }

            if (f->isDiffuse()) { // not possible for mirrors/glass
                Color3f Ld = directLightSampling(scene, sampler, r, f, its);
                if (Ld.isValid()) res += Ld * throughput;
            }

            // BSDF
            BSDFQueryRecord rec(its.shFrame.toLocal(-r.d));
            Color3f F = f->sample(rec, sampler->next2D());

            /*
            When rays from the camera are refracted to different media, the radiance
            they carry is scaled depending on the relative indices of refraction.It's
            worthwhile to apply Russian roulette without including this scaling, since
            it saves us from terminating rays that are actually about to refract out of
            a medium and have the effect of radiance scaling (due to relative IOR) cancelled out.
            */
            eta *= 1.f/rec.eta;
            throughput *= F;
            specular = !f->isDiffuse();
            r = Ray3f(its.p, its.shFrame.toWorld(rec.wo));

            // Russian Roulette
            float pRR = 1.f;
            if (bounces > 3) {
                pRR = std::min(throughput.maxCoeff() * eta * eta, 0.99f);
            }
            throughput /= pRR;

            if (sampler->next1D() >= pRR) break;
            if (!throughput.isValid()) break;
        }

        return res;
    }

    Color3f directLightSampling(const Scene* scene, Sampler* sampler, const Ray3f& ray, const BSDF* f,
        const Intersection& its) const {
        // Uniform-randomly pick a light.
        const std::vector<Mesh*>& lights = scene->getEmitters();
        const Mesh* light = lights[sampler->nextUInt(lights.size())];

        // Sample the light
        EmitterQueryRecord eRec;
        eRec.origin = its.p;
        eRec.origN = its.shFrame.n;
        eRec.light = light;
        Color3f Ld = light->getEmitter()->sample(sampler->next2D(), eRec);
        if (!Ld.isValid()) return Color3f(0.f);
        Ld /= (1.f / lights.size()); // Prob of choosing light

        // Visibility
        Vector3f v = eRec.hit - eRec.origin;
        float d = v.norm();
        v /= d;
        Ray3f shadow(eRec.origin, v, Epsilon, d - Epsilon);
        if (scene->rayIntersect(shadow)) return Color3f(0.f);

        // Direct lighting
        BSDFQueryRecord bRecDirect(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(v), ESolidAngle);
        Color3f Fd = f->eval(bRecDirect);
        Ld *= Fd;

        if (!Ld.isValid()) return Color3f(0.f);

        return Ld;
    }

    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return tfm::format(
            "PathEms"
        );
    }
};

NORI_REGISTER_CLASS(PathEms, "path_ems");
NORI_NAMESPACE_END