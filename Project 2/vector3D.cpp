#include "vector3D.hpp"

Vector3::Vector3() : x(0), y(0), z(0) {}
Vector3::Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

// Addition
Vector3 Vector3::operator+(const Vector3& other) const {
    return Vector3(x + other.x, y + other.y, z + other.z);
}

// Subtraction
Vector3 Vector3::operator-(const Vector3& other) const {
    return Vector3(x - other.x, y - other.y, z - other.z);
}

// Scalar multiplication
Vector3 Vector3::operator*(double scalar) const {
    return Vector3(x * scalar, y * scalar, z * scalar);
}

// Scalar division
Vector3 Vector3::operator/(double scalar) const {
    return Vector3(x / scalar, y / scalar, z / scalar);
}

// Dot product
double Vector3::dot(const Vector3& other) const {
    return x * other.x + y * other.y + z * other.z;
}

// Norm (length)
double Vector3::norm() const {
    return std::sqrt(x * x + y * y + z * z);
}

// Norm squared
double Vector3::norm2() const {
    return x * x + y * y + z * z;
}

// Print
void Vector3::print() const {
    std::cout << "(" << x << ", " << y << ", " << z << ")\n";
}