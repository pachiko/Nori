#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN


// Direct lighting (soft shadows)
class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList& props) {
        m_lightPos = props.getPoint("position");
        // Point ligh has isotropic radiation. so albedo/pi/4pi
        m_lightColor = props.getColor("energy")*INV_FOURPI/M_PI;
    }

    /// Compute the radiance value for a given ray. Just return green here
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        // L V cos(theta) / r^2
        Vector3f delta = m_lightPos - its.p;
        float dist = delta.norm();
        Vector3f v = delta.normalized();

        float cos_theta = v.dot(its.shFrame.n);
        if (cos_theta <= 0) return Color3f(0.0f);

        Ray3f shadow(its.p, v);
        if (scene->rayIntersect(shadow)) return Color3f(0.0f);

        return m_lightColor*cos_theta/dist/dist; // squared-distance fall-off
    }

    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return tfm::format(
            "SimpleIntegrator[\n"
            "  m_lightPos = %s\n"
            "  m_lightColor = %s\n"
            "]",
            m_lightPos,
            m_lightColor
        );
    }
protected:
    Point3f m_lightPos;
    Color3f m_lightColor;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END