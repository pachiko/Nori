/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        Color3f res = m_kd * INV_PI;
        const Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Glancing angles
        if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0) return Color3f(0.f);

        float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        float denom = 4.f * Frame::cosTheta(wh) * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo);
        float G = geometry(bRec);
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);

        res += m_ks * F * G * D / denom;

        return res;
    }

    // Shadowing-Masking function
    // Beckmann-Spizzichino function 
    // 0.5 * [erf(a) - 1 + exp(-a^2)/a/sqrt(PI)]
    // where a = 1/(alpha*tan_theta)
    float geometry(const BSDFQueryRecord& bRec) const {
        const Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        auto lambda = [&](const Vector3f& wv) {
            float chi = wv.dot(wh) / Frame::cosTheta(wv);
            if (chi <= 0.f) return 0.f;

            float absTanTheta = std::abs(Frame::tanTheta(wv));
            float b = 1.f / (m_alpha * absTanTheta);
        
            if (b >= 1.6f) return 1.f;

            return (3.535f * b + 2.181f * b * b) /
                (1.f + 2.276f*b + 2.577f * b * b);
        };

        return lambda(bRec.wo) * lambda(bRec.wi);
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        const Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Glancing angles
        if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0) return 0.f;

        return m_ks * (Warp::squareToBeckmannPdf(wh, m_alpha) / 4.f / wh.dot(bRec.wo))
            + (1.f - m_ks) * (Warp::squareToCosineHemispherePdf(bRec.wo));
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float zeta1 = sinf(sample.x())*10000.f; // https://thebookofshaders.com/10/
        zeta1 -= floor(zeta1); // use modf() alternatively
        bool specular = zeta1 <= m_ks;

        if (specular) { // Specular
            Vector3f wh = Warp::squareToBeckmann(sample, m_alpha); // Normal = halfway
            bRec.wo = 2.f*(wh.dot(bRec.wi))*wh - bRec.wi; // Reflect vector
            bRec.eta = 1.f;
        }
        else { // Diffuse
            bRec.wo = Warp::squareToCosineHemisphere(sample);
            bRec.eta = 1.f;
        }

        Color3f res = eval(bRec);
        if (res.getLuminance() <= 0.f) return Color3f(0.f);
        return res * Frame::cosTheta(bRec.wo) / pdf(bRec);
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
