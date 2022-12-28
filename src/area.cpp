#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN


// Area Light
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
        float prob = pdf(rec) * surfaceRec.pdf;

        if (prob <= 0.f) return Color3f(0.f);
        return res / prob;
    }

    // Radiance
    Color3f eval(EmitterQueryRecord& rec) const {
        Vector3f wi = (rec.origin - rec.hit).normalized();
        Vector3f n = rec.hitN;
        float cos_theta = wi.dot(n);
        return (cos_theta <= 0.f) ? Color3f(0.f) : m_radiance;
    }

    // PDF
    float pdf(EmitterQueryRecord& rec) const {
        Vector3f v = rec.hit - rec.origin;
        float r2 = v.squaredNorm();
        v.normalize();

        float cos_theta_x = v.dot(rec.origN);
        float cos_theta_y = -v.dot(rec.hitN);

        if (cos_theta_x <= 0.f || cos_theta_y <= 0.f) return 0.f;

        return r2 / (cos_theta_x * cos_theta_y);
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