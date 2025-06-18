#pragma once

#include <cmath>
#include <iostream>

class Vector3 {
public:
    double x, y, z;

    Vector3();
    Vector3(double x, double y, double z);

    // Basic operations
    Vector3 operator+(const Vector3& other) const;
    Vector3 operator-(const Vector3& other) const;
    Vector3 operator*(double scalar) const;
    Vector3 operator/(double scalar) const;

    double dot(const Vector3& other) const;
    double norm() const;
    double norm2() const;

    void print() const;
};
