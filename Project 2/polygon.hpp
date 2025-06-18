#pragma once
#include <vector>
#include <iostream>
#include "vector.hpp" 
#include "lbfgs.h"

class MyPolygon {
public:
    std::vector<Vector> vertices;
    
    void print(const std::string& label = "") const {
        if (!label.empty()) std::cout << label << ":\n";
        for (size_t i = 0; i < vertices.size(); ++i) {
            std::cout << "  (" << vertices[i][0] << ", " << vertices[i][1] << ")\n";
        }
        std::cout << "---\n";
    }
};

Vector Intersect(Vector prevVertex, Vector curVertex, Vector edgeStart, Vector edgeEnd);

MyPolygon Sutherland(MyPolygon ClipPolygon, MyPolygon SubjectPolygon);

MyPolygon Sutherland_VPLE(MyPolygon, Vector, Vector);
MyPolygon VPLE(std::vector<Vector>, Vector, MyPolygon);

MyPolygon weighted_Sutherland_VPLE(MyPolygon subjectPolygon, Vector Pi, Vector Pj, double wi, double wj);
MyPolygon weighted_VPLE(
    const std::vector<Vector>& points,
    const std::vector<double>& weights,
    int origin_idx,
    const MyPolygon& Box);

double compute_area(const MyPolygon& poly);
double polygon_integral_squared_distance(const MyPolygon& poly, const Vector& Pi);


struct EvaluateData {
    const std::vector<Vector>* points;
    const MyPolygon* box;
    const std::vector<double>* lambda;
};


lbfgsfloatval_t evaluate(
    void *instance_data,
    const lbfgsfloatval_t *w, 
    lbfgsfloatval_t *g, 
    int n,
    lbfgsfloatval_t step
);

void gallouet_merigot_step(
    std::vector<Vector>& X,
    std::vector<Vector>& v,
    const std::vector<double>& m,
    const MyPolygon& Box,
    const std::vector<double>& lambda,
    double dt,
    double epsilon,
    const Vector& gravity);

MyPolygon clip_with_disk(const MyPolygon& cell, const Vector& center, double radius, int sides = 16);
