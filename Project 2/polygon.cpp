#include "polygon.hpp"
#include "vector.hpp"
#include "vector3D.hpp"
#include "lbfgs.h"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>
#include <sstream>

#include "svg.hpp"

Vector Intersect(Vector prevVertex, Vector curVertex, Vector edgeStart, Vector edgeEnd){
    Vector d1 = curVertex - prevVertex;
    Vector d2 = edgeEnd - edgeStart;

    float crossprod = cross(d1,d2);
    if (crossprod == 0){
        return Vector(999,999);
    }

    Vector delta = edgeStart - prevVertex;
    float t = cross(delta,d2) / crossprod;

    return prevVertex + d1 * t;
}

MyPolygon Sutherland(MyPolygon ClipPolygon, MyPolygon SubjectPolygon){
    if (SubjectPolygon.vertices.size() < 3) {
        return MyPolygon(); // return empty
    }
    for (int i = 0; i < ClipPolygon.vertices.size(); i++){
        MyPolygon outPolygon;
        Vector edgeStart = ClipPolygon.vertices[i];
        Vector edgeEnd = ClipPolygon.vertices[(i + 1) % ClipPolygon.vertices.size()];
        //std::cout << edgeStart[0] << " " << edgeStart[1] << " & " << edgeEnd[0] << edgeEnd[1] << std::endl;
        for (int j = 0; j < SubjectPolygon.vertices.size(); j++){
            Vector curVertex = SubjectPolygon.vertices[j];
            Vector prevVertex = SubjectPolygon.vertices[(j > 0)?(j - 1):(SubjectPolygon.vertices.size()-1)];
            Vector intersection = Intersect(prevVertex, curVertex, edgeStart, edgeEnd);
            if (intersection[0] == 999 && intersection[1] == 999) continue;
            if (cross(edgeEnd - edgeStart, curVertex - edgeStart) >= 0) {
                if (cross(edgeEnd - edgeStart, prevVertex - edgeStart) <= 0) {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if (cross(edgeEnd - edgeStart, prevVertex - edgeStart) >= 0){
                outPolygon.vertices.push_back(intersection);
            }
        }
        SubjectPolygon = outPolygon;
    }
    
    return SubjectPolygon;
}




MyPolygon Sutherland_VPLE(MyPolygon subjectPolygon, Vector Pi, Vector Pj) {
    MyPolygon outPolygon;
    Vector M = (Pi + Pj) * 0.5;
    Vector diff = Pj - Pi;
    Vector normal = Vector(diff[1], -diff[0]);

    Vector edgeStart = M + normal * 100;
    Vector edgeEnd = M - normal*100;
    for (int i = 0; i < subjectPolygon.vertices.size(); i++) {
        Vector curVertex = subjectPolygon.vertices[i];
        Vector prevVertex = subjectPolygon.vertices[(i > 0)?(i - 1):(subjectPolygon.vertices.size()-1)];
        //std::cout << prevVertex[0] << " " << prevVertex[1] << " & " << curVertex[0] << " " << curVertex[1] << std::endl;
        Vector intersection = Intersect(prevVertex, curVertex, edgeStart, edgeEnd);
        if (intersection[0] == 999 && intersection[1] == 999) continue;
        if (cross(edgeEnd - edgeStart, curVertex - edgeStart) >= 0) {
            if (cross(edgeEnd - edgeStart, prevVertex - edgeStart) <= 0) {
                outPolygon.vertices.push_back(intersection);
            }
            outPolygon.vertices.push_back(curVertex);
        }
        else if (cross(edgeEnd - edgeStart, prevVertex - edgeStart) >= 0){
            outPolygon.vertices.push_back(intersection);
        }
    }

    return outPolygon;
}


MyPolygon VPLE(std::vector<Vector> points, Vector origin, MyPolygon Box){
    MyPolygon cell = Box;
    std::sort(points.begin(), points.end(), [&origin](const Vector& a, const Vector& b) {
    return (a - origin).norm2() < (b - origin).norm2();
    });

    double maxDistance = 9999;
    #pragma omp parallel
    for(int i = 0; i < points.size(); i++){
        if (points[i][0] == origin[0] && points[i][1] == origin[1]) continue;

        if ((points[i] - origin).norm() > maxDistance) break;

        cell = Sutherland_VPLE(cell, origin, points[i]);

        double farthest = 0.0;
        for (const Vector& v : cell.vertices) {
            double d = (v - origin).norm();
            if (d > farthest) farthest = d;
        }

        maxDistance = 2.0 * farthest;
    }

    return cell;
}


MyPolygon weighted_Sutherland_VPLE(MyPolygon subjectPolygon, Vector Pi, Vector Pj, double wi, double wj) {
    MyPolygon outPolygon;
    Vector M = (Pi + Pj) * 0.5;
    Vector diff = Pj - Pi;
    if (diff.norm2() < 1e-12) return subjectPolygon;
    Vector M_weighted = M + (wj - wi) * (Pj - Pi) / (2 * diff.norm2());

    Vector normal = Vector(diff[1], -diff[0]);

    Vector edgeStart = M_weighted + normal * 100;
    Vector edgeEnd = M_weighted - normal*100;
    for (int i = 0; i < subjectPolygon.vertices.size(); i++) {
        Vector curVertex = subjectPolygon.vertices[i];
        Vector prevVertex = subjectPolygon.vertices[(i > 0)?(i - 1):(subjectPolygon.vertices.size()-1)];
        //std::cout << prevVertex[0] << " " << prevVertex[1] << " & " << curVertex[0] << " " << curVertex[1] << std::endl;
        Vector intersection = Intersect(prevVertex, curVertex, edgeStart, edgeEnd);
        if (intersection[0] == intersection[1] == 999) {continue;}
        if (cross(edgeEnd - edgeStart, curVertex - edgeStart) >= 0) {
            if (cross(edgeEnd - edgeStart, prevVertex - edgeStart) <= 0) {
                outPolygon.vertices.push_back(intersection);
            }
            outPolygon.vertices.push_back(curVertex);
        }
        else if (cross(edgeEnd - edgeStart, prevVertex - edgeStart) >= 0){
            outPolygon.vertices.push_back(intersection);
        }
    }

    return outPolygon;
}


MyPolygon weighted_VPLE(
    const std::vector<Vector>& points,
    const std::vector<double>& weights,
    int origin_idx,
    const MyPolygon& Box)
{
    MyPolygon cell = Box;
    Vector Pi = points[origin_idx];
    double wi = weights[origin_idx];

    double m = *std::max_element(weights.begin(), weights.end());

    std::vector<int> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { //lambda function to sort points by power distance
        if (i == origin_idx) return false;
        if (j == origin_idx) return true;
        double di = (points[i] - Pi).norm2() - weights[i];
        double dj = (points[j] - Pi).norm2() - weights[j];
        return di < dj;
    });

    double maxDist3D = 9999;

    for (int idx : indices) {
        if (idx == origin_idx) continue;
        Vector Pj = points[idx];
        double wj = weights[idx];

        Vector3 Pi_3D(Pi[0], Pi[1], std::sqrt(std::max(m - wi, 0.0)));
        Vector3 Pj_3D(Pj[0], Pj[1], std::sqrt(std::max(m - wj, 0.0)));

        if ((Pj_3D - Pi_3D).norm() > maxDist3D)
            break;

        cell = weighted_Sutherland_VPLE(cell, Pi, Pj, wi, wj);

        double max_v_dist = 0.0;
        for (const Vector& v : cell.vertices) {
            Vector3 v3D(v[0], v[1], 0.0);
            double dist = (v3D - Pi_3D).norm();
            if (dist > max_v_dist) max_v_dist = dist;
        }

        maxDist3D = 2.0 * max_v_dist;
    }

    return cell;
}

double compute_area(const MyPolygon& poly) {
    double area = 0.0;
    for (int i = 0; i < poly.vertices.size(); ++i) {
        const Vector& p1 = poly.vertices[i];
        const Vector& p2 = poly.vertices[(i + 1) % poly.vertices.size()];
        area += .5 * (p1[0] * p2[1] - p2[0] * p1[1]);
    }
    return std::abs(area) ;
}

double polygon_integral_squared_distance(const MyPolygon& poly, const Vector& Pi) {
    double out = 0.0;
    int n = poly.vertices.size();

    for (int k = 0; k < n; ++k) {
        const Vector& v0 = poly.vertices[k];
        const Vector& v1 = poly.vertices[(k + 1) % n];

        double x0 = v0[0], y0 = v0[1];
        double x1 = v1[0], y1 = v1[1];

        double a = x0 * y1 - x1 * y0;

        double b = x0 * x0 + x0 * x1 + x1 * x1 +  y0 * y0 + y0 * y1 + y1 * y1;

        double c = -4.0 * (Pi[0] * (x0 + x1) + Pi[1] * (y0 + y1)) - 6 * Pi.norm2();

        out +=  a * (b + c);
    }

    return std::abs(out) / 12.0;
}

lbfgsfloatval_t evaluate(
    void *instance_data,
    const lbfgsfloatval_t *w,
    lbfgsfloatval_t *g,
    int n,
    lbfgsfloatval_t step
) {
    auto* data = static_cast<EvaluateData*>(instance_data);
    const std::vector<Vector>& points = *(data->points);
    const MyPolygon& Box = *(data->box);
    const std::vector<double>& lambda = *(data->lambda);

    double total_energy = 0.0;

    for (int i = 0; i < n; ++i) {
        MyPolygon cell = weighted_VPLE(points, std::vector<double>(w, w + n), i, Box);
        if (cell.vertices.size() < 3) {
            g[i] = 0.0;
            continue;}
        double area = compute_area(cell);
        double integral = polygon_integral_squared_distance(cell, points[i]);

        total_energy += integral + lambda[i] * w[i] - w[i] * area;
        g[i] = (area - lambda[i]);
    }

    return -total_energy;
}


void gallouet_merigot_step(
    std::vector<Vector>& X,
    std::vector<Vector>& v,
    const std::vector<double>& m,
    const MyPolygon& Box,
    const std::vector<double>& lambda,
    double dt,
    double epsilon,
    const Vector& gravity)
{
    const int N = X.size();
    std::vector<double> weights(N, 0.0);

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.max_iterations = 500;

    EvaluateData data;
    data.points = &X;
    data.box = &Box;
    data.lambda = &lambda;

    lbfgs(N, weights.data(), nullptr, evaluate, nullptr, &data, &param);

    std::vector<Vector> new_X(N), new_v(N);

    for (int i = 0; i < N; ++i) {
        MyPolygon cell = weighted_VPLE(X, weights, i, Box);
        if (cell.vertices.size() < 3) {
            new_X[i] = X[i];
            new_v[i] = v[i];
            continue;
        }

        Vector centroid(0, 0);
        double total_area = 0.0;
        for (int j = 0; j < cell.vertices.size(); ++j) {
            
            const Vector& p1 = cell.vertices[j];
            const Vector& p2 = cell.vertices[(j + 1) % cell.vertices.size()];
            double a = p1[0] * p2[1] - p2[0] * p1[1];
            Vector tri_center = (p1 + p2) * (1.0 / 3.0);
            centroid = centroid + tri_center * a;
            total_area += a;
        }
        centroid = centroid * (1.0 / total_area);

        Vector F_spring = (centroid - X[i]) * (1.0 / (epsilon * epsilon));
        Vector F_total = F_spring + gravity * m[i];

        new_v[i] = v[i] + F_total * (dt / m[i]);
        new_X[i] = X[i] + new_v[i] * dt;
    }

    X = std::move(new_X);
    v = std::move(new_v);
}

MyPolygon clip_with_disk(const MyPolygon& cell, const Vector& center, double radius, int sides) {
    MyPolygon disk;
    for (int i = 0; i < sides; ++i) {
        double angle = 2.0 * M_PI * i / sides;
        disk.vertices.push_back(center + Vector(std::cos(angle), std::sin(angle)) * radius);
    }
    return Sutherland(disk, cell);
}