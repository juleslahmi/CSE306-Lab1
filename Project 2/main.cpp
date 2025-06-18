#include "svg.hpp"
#include "vector.hpp"
#include "polygon.hpp"
#include "lbfgs.h"
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <fstream>

int main(){
    /*
    //---------------------------------------------------------------------------------
    //Test of SutherLand

    std::vector<MyPolygon> polygons;

    MyPolygon clipPolygon;
    clipPolygon.vertices.push_back(Vector(0.2, 0.2));
    clipPolygon.vertices.push_back(Vector(0.8, 0.2));
    clipPolygon.vertices.push_back(Vector(0.6, 0.8));
    clipPolygon.vertices.push_back(Vector(0.4, 0.8));

    MyPolygon subjectPolygon;
    subjectPolygon.vertices.push_back(Vector(0.3, 0.3));
    subjectPolygon.vertices.push_back(Vector(0.9, 0.3));
    subjectPolygon.vertices.push_back(Vector(0.9, 0.9));
    subjectPolygon.vertices.push_back(Vector(0.3, 0.9));

    polygons.push_back(clipPolygon);
    polygons.push_back(subjectPolygon);


    save_svg(polygons, "polygons.svg");
    MyPolygon output = Sutherland(clipPolygon, subjectPolygon);
    for (int i = 0; i < output.vertices.size(); i++){
        //std::cout << output.vertices[i][0] << " " << output.vertices[i][1] << std::endl;
    }
    std::vector<MyPolygon> polygons_2;
    polygons_2.push_back(output);
    save_svg(polygons_2, "polygons_clipped.svg");
    */

    //---------------------------------------------------------------------------------
    //Test of VPLE
    /*
    std::vector<Vector> points;
    std::vector<double> target_masses;

    std::default_random_engine gen;
    std::uniform_real_distribution<double> dis(0, 1); 

    for (int i = 0; i < 2000; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        points.emplace_back(x, y);

    }

    std::vector<MyPolygon> diagram;
    MyPolygon box;
    box.vertices = {
        Vector(0.0, 0.0), Vector(1.0, 0.0),
        Vector(1.0, 1.0), Vector(0.0, 1.0)
    };

    for (int i = 0; i < points.size(); ++i) {
        diagram.push_back(VPLE(points, points[i], box));
    }
    save_svg(diagram, "VPLE.svg");
    */

    //---------------------------------------------------------------------------------
    //Test of Weighted VPLE

    /*
    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<double> target_masses;

    std::default_random_engine gen;
    std::uniform_real_distribution<double> dis(0, 1); 

    for (int i = 0; i < 2000; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        points.emplace_back(x, y);
        
        Vector p(x, y);
        double d2 = (p - Vector(0.5, 0.5)).norm2();
        weights.push_back(2 * std::exp(-d2 / 0.02));
    }

    std::vector<MyPolygon> diagram;
    MyPolygon box;
    box.vertices = {
        Vector(0.0, 0.0), Vector(1.0, 0.0),
        Vector(1.0, 1.0), Vector(0.0, 1.0)
    };

    for (int i = 0; i < points.size(); ++i) {
        diagram.push_back(weighted_VPLE(points, weights, i, box));
    }
    save_svg(diagram, "Gaussian_Spread.svg");
    */
        
    //---------------------------------------------------------------------------------
    //Test of Gallouet Algorithm

    /*
    MyPolygon Box;
    Box.vertices = {
        Vector(0, 0),
        Vector(1, 0),
        Vector(1, 1),
        Vector(0, 1)
    };

    std::default_random_engine gen;
    std::uniform_real_distribution<double> dis(0, 1); 
    std::vector<Vector> X;
    for (int i = 0; i < 100; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        X.emplace_back(x, y);
    }
    std::vector<Vector> v(X.size(), Vector(0, 0));
    std::vector<double> m(X.size(), 1.0);
    std::vector<double> lambda(X.size(), 1.0 / X.size());
    Vector gravity(0, 0);
    double dt = 0.01;
    double epsilon = 0.05;

    std::vector<double> weights(X.size(), 0.0);

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.max_iterations = 100;

    EvaluateData data;
    data.points = &X;
    data.box = &Box;
    data.lambda = &lambda;

    lbfgs(X.size(), weights.data(), nullptr, evaluate, nullptr, &data, &param);

    std::vector<MyPolygon> test;
    for (size_t i = 0; i < X.size(); ++i) {
        MyPolygon cell = weighted_VPLE(X, weights, i, Box);
        test.push_back(cell);
    }
    save_svg(test, "Gallouet_Before.svg");

    for (int i = 0; i < 100; i++){
    gallouet_merigot_step(X, v, m, Box, lambda, dt, epsilon, gravity);
    }
    for (size_t i = 0; i < X.size(); ++i) {
        std::cout << "Particle " << i << " new position: " << X[i][0] << ", " << X[i][1]
                << " | new velocity: " << v[i][0] << ", " << v[i][1] << std::endl;
    }

    std::vector<double> weights_after(X.size(), 0.0);
    data.points = &X;
    lbfgs(X.size(), weights_after.data(), nullptr, evaluate, nullptr, &data, &param);

    std::vector<MyPolygon> test1;
    for (size_t i = 0; i < X.size(); ++i) {
        MyPolygon cell = weighted_VPLE(X, weights_after, i, Box);
        test1.push_back(cell);
    }
    save_svg(test1, "Gallouet_After.svg");
    return 0;
    */
    //---------------------------------------------------------------------------------
    //Test of Fluid Mechanics (does not work!)

    MyPolygon Box;
    Box.vertices = {
        Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1)
    };

    run_fluid_simulation(
        25,          // number of fluid particles
        100,         // number of frames
        0.01,        // time step
        0.05,        // spring stiffness
        0.6,         // total fluid volume
        Box,         // simulation domain
        Vector(0, -9.81)  // gravity
    );

    return 0;
}

