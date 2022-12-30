#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>

NORI_NAMESPACE_BEGIN


// Whitted
class Whitted : public Integrator {
public:
    Whitted(const PropertyList& props) { }

    /// Compute the radiance value for a given ray. Just return green here
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection first;
        if (!scene->rayIntersect(ray, first)) return Color3f(0.f); // No hit

        const BSDF* f = first.mesh->getBSDF();

        if (f->isDiffuse()) {
            // Hit a light?
            Color3f Le(0.f);
            if (first.mesh->isEmitter()) {
                EmitterQueryRecord rec;
                rec.origin = ray.o;
                rec.hit = first.p;
                rec.hitN = first.shFrame.n;
                rec.light = first.mesh;
                Le = first.mesh->getEmitter()->eval(rec);
            }

            // Uniform-randomly pick a light.
            const std::vector<Mesh*>& lights = scene->getEmitters();
            const Mesh* light = lights[sampler->nextUInt(lights.size())];

            // Sample the light
            EmitterQueryRecord rec;
            rec.origin = first.p;
            rec.origN = first.shFrame.n;
            rec.light = light;
            Color3f Lr = light->getEmitter()->sample(sampler->next2D(), rec);
            if (Lr.getLuminance() <= 0.f) return Le;
            Lr /= (1.f / lights.size()); // Prob of choosing light

            // Visibility
            Vector3f v = rec.hit - rec.origin;
            float d = v.norm();
            v /= d;
            Ray3f shadow(rec.origin, v, Epsilon, d - Epsilon);
            if (scene->rayIntersect(shadow)) return Le;

            // BSDF
            BSDFQueryRecord bRec(first.shFrame.toLocal(-ray.d), first.shFrame.toLocal(v), ESolidAngle);
            Color3f F = f->eval(bRec);
            Lr *= F;

            return Le + Lr;
        } else if (sampler->next1D() < 0.95f) {
            BSDFQueryRecord bRec(first.shFrame.toLocal(-ray.d));
            Color3f res = f->sample(bRec, sampler->next2D());

            Ray3f next(first.p, first.shFrame.toWorld(bRec.wo));
            return res*Li(scene, sampler, next)/0.95f;
        }

        return Color3f(0.f);
    }

    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return tfm::format(
            "Whitted"
        );
    }
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END