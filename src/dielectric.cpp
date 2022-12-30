/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        bRec.measure = EDiscrete;

        float cosThetaI = Frame::cosTheta(bRec.wi);
        float Fr = fresnel(cosThetaI, m_extIOR, m_intIOR);

        if (sample.x() < Fr) {
            // reflect
            bRec.eta = 1.f;

            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                cosThetaI
            );

            return Color3f(1.f);
        } else { 
            // refract
            float etaI = m_extIOR, etaT = m_intIOR;

            /* Swap the indices of refraction if the interaction starts
            at the inside of the object */
            if (cosThetaI < 0.0f) {
                std::swap(etaI, etaT);
                // No need to negate cosThetaI since we square it here.
            }

            /* Using Snell's law, calculate the squared sine of the
            angle between the normal and the transmitted ray */
            float eta = etaI / etaT,
                sinThetaTSqr = eta * eta * (1 - cosThetaI * cosThetaI);
            float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

            // what goes up, goes down; vice-versa.
            if (cosThetaI > 0.f) cosThetaT = -cosThetaT;

            bRec.eta = eta;

            // vector component parallel to normal is cos, so perpendicular is sin
            // use Snell's Law: eta_i sin_i = eta_t sin_t
            bRec.wo = Vector3f(
                -bRec.wi.x()*eta,
                -bRec.wi.y()*eta,
                cosThetaT
            ).normalized();

            return eta*eta*Color3f(1.f);
        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
