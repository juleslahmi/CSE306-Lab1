#include "vector.hpp"

Vector::Vector(double x, double y) {
    data[0] = x;
    data[1] = y;
}

double Vector::norm2() const {
    return data[0] * data[0] + data[1] * data[1];
}

double Vector::norm() const {
    return std::sqrt(norm2());
}

void Vector::normalize() {
    double n = norm();
    if (n > 0) {
        data[0] /= n;
        data[1] /= n;
    }
}

double Vector::operator[](int i) const {
    return data[i];
}
double& Vector::operator[](int i) {
    return data[i];
}

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a * b[0], a * b[1]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0] * b, a[1] * b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1];
}

double cross(const Vector& a, const Vector& b){
    return a[0] * b[1] - a[1] * b[0];
}
