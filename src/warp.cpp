/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

// Inversion Sampling
// 1. Make sure the pdf is really a pdf (does it integrate to 1?)
// 1a. If not, integrate c.f(x) = 1. Solve for c, and p(x) = c.f(x)
// 2. Compute CDF: P(x) by integrating p(x') for 0 < x' < x
// 3. Compute inverse CDF P-1(x) by rearranging P(y) = x
// 4. Compute x = P-1(zeta) (Mapped RV), 0 < zeta < 1 (Uniform RV)
// 5. Use x to compute p(x) (PDF)


Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    // 1. BREAK p(t) = 1 - |t| into -1 < t < 0 and 0 < t < 1
    // 2. integrate p(t) = 1 + t for -1 < t < 0, to obtain P(0) = 1/2. P(-1) = 0.
    // 3. integrate p(t') = 1 + t for -1 < t' < t, to obtain P(t) = t + t * t / 2.f + 0.5f
    // 3a. verify P(0) = 1/2, and P(-1) = 0
    // 4. integrate p(t') = 1 - t for 0 < t' < t, to obtain P(t) = t - t * t / 2.f
    // 4a. ADD P(0) into P(t) to obtain P(t) = t - t * t / 2.f + 0.5f
    // 4b. verify P(0) = 1/2, and P(1) = 1

    /*
    auto tentCDF = [](float t){
        if (t >= -1 && t <= 0) return t + t * t / 2.f + 0.5f;
        else if (t >= 0 && t <= 1) return t - t * t / 2.f + 0.5f;
        else return 0.f;
    };
    */

    // Given zeta and using CDF, solve for t. Requires solving a quadratic equation.
    auto tentInv = [](float z) {
        if (z >= 0.f && z <= 0.5f) {
            return -1.f + sqrtf(2.f * z);
        }
        else if (z >= 0.5f && z <= 1.f) {
            return 1.f - sqrtf(2.f - 2.f * z);
        }
        else {
            return -1.f;
        }
    };

    float tx = tentInv(sample.x());
    float ty = tentInv(sample.y());

    return Point2f(tx, ty);
}

float Warp::squareToTentPdf(const Point2f &p) {
    auto tentPDF = [](float t) {
        if (t >= -1.f && t <= 0.f) {
            return 1.f + t;
        }
        else if (t >= 0.f && t <= 1.f) {
            return 1.f - t;
        }
        else {
            return 0.f;
        }
    };

    return tentPDF(p.x()) * tentPDF(p.y());
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    // 1. Area patch of unit disk is dA = r * theta d_theta d_r
    //    Integrate to obtain A = r^2/2 * 2PI = PI. Invert to obtain p(x, y) = 1 / PI
    // 2. transform to p(r, theta), given x = r*cos(theta), y = r*sin(theta)
    // 3. p(y) = |dy/dx|^(-1)*p(x)
    // 4. we are drawing samples from p(r, theta), so p(x, y) = 1/|del(x, y)/del(r, theta)|*p(r, theta)
    //    |del(x, y)/del(r, theta)| = | del_x/del_r, del_x/del_theta \ del_y/del_r, del_y/del_theta |
    //    |(cos(theta) * r cos(theta)) - (-r sin(theta) * sin(theta))| = r(cos^2 + sin^2) = r
    //    1/PI = 1/r*p(r, theta); p(r, theta) = r/PI
    // 5. now that we have the joint PDF, integrate one of the variables to get the marginal
    //    p(theta) = int(1->0)[r/pi]dr = [r^2/2pi](1->0) = 1/2pi  OR
    //    p(r) = int(2pi->0)[r/pi]dtheta = [r*theta/pi](2pi->0) = 2r
    // 6. Use Bayes' Theorem to get the conditional p(x|y) = p(x, y)/p(y)
    //    p(r|theta) = 2r, p(theta|r) = 1/2pi; coincidence?
    // 7. Integrate to get CDF, then invert
    //    P(r) = r^2, zeta = r^2, r = sqrt(zeta)
    //    P(theta|r) = theta/2pi, zeta = theta/2pi, theta = 2pi*zeta

    float r = sqrtf(sample.x());
    float theta = 2 * M_PI * sample.y();
    return Point2f(r*cos(theta), r*sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    // Sampling proportional to area of unit disk, so p(x, y) = 1/PI
    if (p.norm() > 1.f) return 0;
    return INV_PI;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    // Area patch for sphere is dA = r^2 * sin(theta) d_theta d_phi
    // Alternatively, use the determinant of Jacobian (3x3)
    // Since dA = dw/r^2, dA = dw = sin(theta) d_theta d_phi 
    // Integrate to obtain A = r^2[-cos(PI) + cos(0)][2PI] = 4PI.
    // Invert to obtain p(w) = 1/4PI
    // p(theta, phi) d_theta d_phi = p(w) dw
    // p(theta, phi) = sin(theta) / 4PI

    // Integrate phi or theta to get marginal:
    // p(theta) = int(2pi->0)[sin(theta)/4pi]d_phi = sin(theta)/2
    // p(phi) = int(pi->0)[sin(theta)/4pi]d_theta = [-cos(theta)/4pi](pi->0) = 1/2pi

    // Bayes law for conditional.
    // p(phi|theta) = sin(theta)/4pi / (sin(theta)/2) = 1/2pi
    // p(theta|phi) = sin(theta)/4pi / (1/2pi) = sin(theta)/2
    // coincidence? I think not.

    // Integrate for CDF:
    // P(theta) = int(theta->0)[sin(theta')/2]d_theta' = [-cos(theta')/2](theta->0) 
    // = [-cos(theta) --cos(0)]/2 = [1 - cos(theta)]/2
    // 1 - cos(theta) = 2*zeta, theta = acos(1 - 2*zeta)
    // P(phi) = int(phi->0)[1/2pi]d_phi' = phi/2pi
    // phi/2pi = zeta, phi = 2pi*zeta

    float theta = acos(1.f - 2.f * sample.x());
    float phi = sample.y() * 2 * M_PI;
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return INV_FOURPI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    // Area patch for sphere is dA = r^2 * sin(theta) d_theta d_phi
    // Alternatively, use the determinant of Jacobian (3x3)
    // Since dA = dw/r^2, dA = dw = sin(theta) d_theta d_phi 
    // Integrate to obtain A = r^2[-cos(pi/2) + cos(0)][2PI] = 2PI.
    // Invert to obtain p(w) = 1/2PI
    // p(theta, phi) d_theta d_phi = p(w) dw {int any pdf = 1}
    // p(theta, phi) = sin(theta) / 2PI

    // Integrate phi or theta to get marginal:
    // p(theta) = int(2pi->0)[sin(theta)/2pi]d_phi = sin(theta)
    // p(phi) = int(pi/2->0)[sin(theta)/2pi]d_theta = [-cos(theta)/2pi](pi/2->0) = 1/2pi

    // Bayes law for conditional.
    // p(phi|theta) = sin(theta)/2pi / sin(theta) = 1/2pi
    // p(theta|phi) = sin(theta)/2pi / (1/2pi) = sin(theta)

    // Integrate for CDF:
    // P(theta) = int(theta->0)[sin(theta')]d_theta' = [-cos(theta')](theta->0) 
    // = [-cos(theta) --cos(0)] = 1 - cos(theta)
    // 1 - cos(theta) = zeta, theta = acos(1 - zeta)
    // P(phi) = int(phi->0)[1/2pi]d_phi' = phi/2pi
    // phi/2pi = zeta, phi = 2pi*zeta

    float theta = acos(1.f - sample.x());
    float phi = sample.y() * 2 * M_PI;
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    float cos_theta = Frame::cosTheta(v);
    if (cos_theta < 0.f) return 0.f;
    return INV_TWOPI;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    // Make sure the pdf is really a pdf (does it integrate to 1 ?)
    // p(w) prop. cos(theta). p(w) = c*cos(theta)
    // int[p(w)]dw = int(pi/2->0, 2pi->0)[c*cos(theta)*sin(theta)] d_theta d_phi 
    // int(pi/2->0, 2pi->0)[c/2*sin(2*theta)] d_theta d_phi
    // c/2*2pi[-cos(2*theta)/2](pi/2->0) = 1
    // c/2*pi[-cos(2*pi/2) - -cos(0)] = 1
    // c/2*pi[-(-1) + 1] = 1; c = 1/pi
    // p(w) = cos(theta)/pi

    // p(theta, phi) = sin(theta)*cos(theta)/pi = sin(2*theta)/2pi

    // Obtain marginals
    // p(theta) = int(2pi->0)[sin(2*theta)/2pi] d_phi = sin(2*theta)
    // p(phi) = int(pi/2->0)[sin(2*theta)/2pi] d_theta = 1/2pi[-cos(2*theta)/2](pi/2->0)
    // = 1/2pi[-(-1)/2 - (-1)/2] = 1/2pi

    // Obtain conditionals
    // p(theta|phi) = p(theta, phi)/p(phi) = sin(2*theta)/2pi/(1/2pi) = sin(2*theta)
    // p(phi|theta) = p(theta, phi)/p(theta) = sin(2*theta)/2pi/sin(2*theta) = 1/2pi
    
    // Integrate to get CDF:
    // P(theta) = int(theta->0)[sin(2*theta')]d_theta' = [-cos(2*theta')/2](theta->0)
    // = [-cos(2*theta)/2 + cos(0)/2] = (1-cos(2*theta))/2
    // theta = acos(1 - 2*zeta)
    // P(phi) = int(phi->0)[1/2pi]d_phi' = [phi'/2pi](phi->0) = phi/2pi
    // phi = 2pi*zeta

    float theta = acos(1.f - 2.f*sample.x())/2.f;
    float phi = sample.y() * 2 * M_PI;
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float cos_theta = Frame::cosTheta(v);
    if (cos_theta < 0.f) return 0.f;
    return cos_theta * INV_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    // Let x = cos(theta),
    // tan^2(theta) = sin^2(theta)/cos^2(theta) = (1 - cos^2(theta))/cos^2(theta) = (1 - x^2)/x^2
    // d_x = -sin(theta) d_theta

    // Longitudinal: 2*exp(-tan^2(theta)/a^2)*sin(theta)/(a^2*cos^3(theta))

    // Integration by substitution for exponential function
    // int[f'(x)*exp(f(x))]d_x = exp(f(x)) + c

    // Substitution Mapping
    // Let f(x) = -tan^2(theta)/a^2 = (x^2 - 1)/(a^2 * x^2)
    // Use Quotient Rule: f'(x) = 2/(a^2 * x^3)

    // Integrate to obtain CDF
    // int(theta->0)[2*exp(-tan^2(theta)/a^2)*sin(theta)/(a^2*cos^3(theta))]d_theta'
    // = int(x->0)[-f'(x)exp(f(x))](d_x) = [-exp(f(x))](x->0) = [-exp(-tan^2(theta')/a^2)](theta->0)
    // = 1 - exp(-tan^2(theta)/a^2)

    // 1 - exp(-tan^2(theta)/a^2) = zeta
    // log(1 - zeta) = -tan^2(theta)/a^2
    // -a^2*log(1 - zeta) = tan^2(theta)
    // theta = atan(sqrt(-a^2*log(1 - zeta)))

    float alpha_sqr = alpha * alpha;
    float logZeta = log(1.f - sample.x());
    if (isinf(logZeta)) logZeta = 0.f;
    float tan2theta = -alpha_sqr*logZeta;

    // tan^2 = (1 - cos^2)/cos^2 ; tan^2*cos^2 = 1 - cos^2 ; (1 + tan^2)cos^2 = 1
    // cos = sqrt(1/(1 + tan^2))
    float cosTheta = 1 / std::sqrt(1 + tan2theta);
    float sinTheta = sqrtf(std::max(0.f, 1.f - cosTheta * cosTheta)); // sin^2 + cos^2 = 1
    float phi = sample.y() * 2 * M_PI;
    float x = sinTheta * cos(phi);
    float y = sinTheta * sin(phi);
    float z = cosTheta;
    return Vector3f(x, y, z);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float cos_theta = Frame::cosTheta(m);
    if (cos_theta <= 0.f) return 0.f; // Prevent divide by zero

    // D(theta, phi)
    float alpha_sqr = alpha * alpha;
    float cos_theta_cube = cos_theta * cos_theta * cos_theta;

    float tan_theta_sqr = Frame::tanTheta(m);
    tan_theta_sqr = tan_theta_sqr * tan_theta_sqr;

    return INV_TWOPI * 2.f * exp(-tan_theta_sqr/alpha_sqr) / (alpha_sqr*cos_theta_cube);
}

NORI_NAMESPACE_END
