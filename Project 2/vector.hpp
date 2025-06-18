#pragma once
#include <cmath>

class Vector {
public:
    explicit Vector(double x = 0, double y = 0);

    double norm2() const;
    double norm() const;
    void normalize();

    double operator[](int i) const;
    double& operator[](int i);

    double data[2];
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const double a, const Vector& b);
Vector operator*(const Vector& a, const double b);
Vector operator*(const Vector& a, const Vector& b);
Vector operator/(const Vector& a, const double b);
double dot(const Vector& a, const Vector& b);
double cross(const Vector& a, const Vector& b);
