#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>


NORI_NAMESPACE_BEGIN


// Ambient Occlusion
class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList& props) {
        srand(static_cast <unsigned> (time(0)));
    }

    /// Compute the radiance value for a given ray. Just return green here
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)) return Color3f(0.0f);
        
        Point2f sample = sampler->next2D();
        Vector3f wi = Warp::squareToCosineHemisphere(sample);
        Frame f(its.shFrame.n);
        wi = f.toWorld(wi);

        Ray3f shadow(its.p, wi);
        if (scene->rayIntersect(shadow)) return Color3f(0.f);

        // NOT rho/PI*cos_theta, since we need to divide by PDF of cosineHemisphere (cos_theta/PI)
        // which equals rho in the end (NB: albedo for AO is 1)
        return Color3f(1.f);
    }

    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return tfm::format(
            "AOIntegrator"
        );
    }
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END