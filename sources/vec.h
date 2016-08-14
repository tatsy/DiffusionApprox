#ifndef _VEC_H_
#define _VEC_H_

#include <cmath>

class Vec;
Vec operator/(const Vec& v, double d);

class Vec {
public:
    // Public methods
    explicit Vec(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
        : x(x_), y(y_), z(z_) {
    }

    ~Vec() {}

    Vec(const Vec& v)
        : x(v.x), y(v.y), z(v.z) {
    }

    Vec& operator=(const Vec& v) {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
        return *this;
    }

    Vec& operator+=(const Vec& v) {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return *this;
    }

    Vec& operator-=(const Vec& v) {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return *this;
    }

    Vec& operator*=(double s) {
        this->x *= s;
        this->y *= s;
        this->z *= s;
        return *this;
    }

    Vec& operator/=(double d) {
        double s = 1.0 / d;
        return this->operator*=(s);
    }

    Vec operator-() const {
        return Vec(-x, -y, -z);
    }

    double dot(const Vec& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec cross(const Vec& v) const {
        double nx = y * v.z - z * v.y;
        double ny = z * v.x - x * v.z;
        double nz = x * v.y - y * v.x;
        return Vec(nx, ny, nz);
    }

    double norm() const {
        return std::sqrt(this->dot(*this));
    }

    Vec normalized() const {
        return (*this) / this->norm();
    }

    // Public parameters
    double x, y, z;
};

Vec operator+(const Vec& v1, const Vec& v2) {
    Vec ret = v1;
    ret += v2;
    return ret;
}

Vec operator-(const Vec& v1, const Vec& v2) {
    Vec ret = v1;
    ret -= v2;
    return ret;
}

Vec operator*(const Vec& v, double s) {
    Vec ret = v;
    ret *= s;
    return ret;
}

Vec operator*(double s, const Vec& v) {
    Vec ret = v;
    ret *= s;
    return ret;
}

Vec operator/(const Vec& v, double d) {
    Vec ret = v;
    ret /= d;
    return ret;
}


#endif  // _VEC_H_
