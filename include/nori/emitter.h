/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

// See BSDFQueryRecord
struct EmitterQueryRecord {
    Point3f origin;
    Normal3f origN;

    Point3f hit;
    Normal3f hitN;

    const Mesh* light;
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:
    // Sample (Color will be divided by the pdf already)
    virtual Color3f sample(const Point2f& sample, EmitterQueryRecord& rec) const = 0;

    // Radiance
    virtual Color3f eval(EmitterQueryRecord& rec) const = 0;

    // PDF
    virtual float pdf(EmitterQueryRecord& rec) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
};

NORI_NAMESPACE_END
