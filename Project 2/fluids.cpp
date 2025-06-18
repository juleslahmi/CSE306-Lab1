#include <vector>
#include "vector.hpp"
#include "polygon.hpp"
#include "vector3D.hpp"
#include "lbfgs.h"
#include "svg.hpp"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>
#include <sstream>

void gallouet_merigot_step_with_air(
    std::vector<Vector>& X,
    std::vector<Vector>& v,
    const std::vector<double>& m,
    const MyPolygon& Box,
    const std::vector<double>& lambda,
    double dt,
    double epsilon,
    const Vector& gravity,
    std::vector<double>& weights
) {
    const int N = X.size();

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.max_iterations = 100;

    EvaluateData data;
    data.points = &X;
    data.box = &Box;
    data.lambda = &lambda;

    lbfgs(N + 1, weights.data(), nullptr, evaluate_fluid_air, nullptr, &data, &param);

    double air_weight = weights[N];
    std::vector<Vector> new_X(N), new_v(N);

    for (int i = 0; i < N; ++i) {
        MyPolygon cell = fluid_air_VPLE(X, weights, air_weight, i, Box);
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
            double a = 0.5 * (p1[0] * p2[1] - p2[0] * p1[1]);
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


MyPolygon fluid_air_VPLE(
    const std::vector<Vector>& points,
    const std::vector<double>& weights,
    double air_weight,
    int idx,
    const MyPolygon& Box)
{
    MyPolygon cell = Box;
    Vector Pi = points[idx];
    double wi = weights[idx];

    double m = *std::max_element(weights.begin(), weights.end());

    std::vector<int> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        if (i == idx) return false;
        if (j == idx) return true;
        double di = (points[i] - Pi).norm2() - weights[i];
        double dj = (points[j] - Pi).norm2() - weights[j];
        return di < dj;
    });

    double maxDist3D = 9999;
    for (int j : indices) {
        if (j == idx) continue;
        Vector Pj = points[j];
        double wj = weights[j];

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

    double radius = std::sqrt(std::max(wi - air_weight, 0.0));
    cell = clip_with_disk(cell, Pi, radius);

    return cell;
}

lbfgsfloatval_t evaluate_fluid_air(
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

    int N = points.size();
    const double air_weight = w[N]; 
    double total_energy = 0.0;
    double estimated_air_volume = 0.0;

    for (int i = 0; i < N; ++i) {
        MyPolygon cell = fluid_air_VPLE(points, std::vector<double>(w, w + N), air_weight, i, Box);
        if (cell.vertices.size() < 3) {
            g[i] = 0.0;
            continue;
        }

        double area = compute_area(cell);
        double integral = polygon_integral_squared_distance(cell, points[i]);

        total_energy += integral + lambda[i] * w[i] - w[i] * area;
        g[i] = lambda[i] - area;
        estimated_air_volume += area;
    }

    double desired_air_volume = 1.0 - std::accumulate(lambda.begin(), lambda.end(), 0.0);
    double air_term = air_weight * (desired_air_volume - estimated_air_volume);
    total_energy += air_term;

    g[N] = desired_air_volume - estimated_air_volume;

    return -total_energy;
}

void run_fluid_simulation(
    int num_particles,
    int num_steps,
    double dt,
    double epsilon,
    double total_fluid_volume,
    const MyPolygon& domain,
    const Vector& gravity
) {
    std::vector<Vector> X;
    std::vector<Vector> v;
    std::vector<double> m;
    std::vector<double> lambda;
    std::vector<double> weights(num_particles + 1, 0.1);   
    weights[num_particles] = 0.01;         
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dis(0.1, 0.9);
    for (int i = 0; i < num_particles; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        X.emplace_back(x, y);
        v.emplace_back(0, 0);
        m.push_back(1.0);
    }

    lambda.assign(num_particles, total_fluid_volume / num_particles);

    for (int step = 0; step < num_steps; ++step) {

        std::vector<double> weights(num_particles + 1, 0.0);
        EvaluateData data{&X, &domain, &lambda};

        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.max_iterations = 100;

        lbfgs(num_particles + 1, weights.data(), nullptr, evaluate_fluid_air, nullptr, &data, &param);

        double air_weight = weights[num_particles];
        std::vector<MyPolygon> cells;
        for (int i = 0; i < num_particles; ++i) {
            MyPolygon cell = fluid_air_VPLE(X, weights, air_weight, i, domain);
            if (cell.vertices.size() >= 3)
                cells.push_back(cell);
        }

        std::stringstream filename;
        filename << "frame_" << step << ".svg";
        save_svg(cells, filename.str());

        gallouet_merigot_step_with_air(X, v, m, domain, lambda, dt, epsilon, gravity, weights);

    }

    std::cout << "Simulation complete." << std::endl;
}
