#ifndef _RAY_H_
#define _RAY_H_

#include "vec.h"

class Ray {
public:
    // Public methods
    Ray() : org(), dir() {}
    Ray(const Vec& o, const Vec& d) : org(o), dir(d) {}
    Ray(const Ray& r) : org(r.org), dir(r.dir) {}

    ~Ray() {}

    Ray& operator=(const Ray& r) {
        this->org = r.org;
        this->dir = r.dir;
        return *this;
    }

    // Public parameters
    Vec org, dir;
};

#endif  // _RAY_H_
