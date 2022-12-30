#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN


// Diffuse Area Light
class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList& props) {
        m_radiance = props.getColor("radiance");
    }

    // Sample (Color will be divided by the pdf already)
    Color3f sample(const Point2f& sample, EmitterQueryRecord& rec) const {
        const Mesh* light = rec.light;

        MeshSurfaceQueryRecord surfaceRec;
        light->sampleSurface(sample, surfaceRec);

        rec.hit = surfaceRec.p;
        rec.hitN = surfaceRec.n;

        Color3f res = eval(rec);
        res *= geometry(rec);
        return res / pdf(rec);
    }

    // Radiance, see DiffuseAreaLight::L()
    Color3f eval(EmitterQueryRecord& rec) const {
        Vector3f v = rec.origin - rec.hit;
        Vector3f w = v.normalized();
        float cos_theta_y = w.dot(rec.hitN);
        return (cos_theta_y > 0.f) ? m_radiance : Color3f(0.f);
    }

    // PDF
    float pdf(EmitterQueryRecord& rec) const {
        return 1.f/rec.light->surfaceArea();
    }

    // Geometry term w/o visibility
    float geometry(EmitterQueryRecord& rec) const {
        Vector3f v = rec.origin - rec.hit;
        Vector3f w = v.normalized();
        float cos_theta_y = w.dot(rec.hitN);
        float cos_theta_x = -w.dot(rec.origN);
        return std::abs(cos_theta_x) * std::abs(cos_theta_y) / v.squaredNorm();
    }

    /// Return a brief string summary of the instance (for debugging purposes)
    std::string toString() const {
        return tfm::format(
            "AreaLight[\n"
            "  m_radiance = %s\n"
            "]",
            m_radiance
        );
    }

protected:
    Color3f m_radiance;

};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END