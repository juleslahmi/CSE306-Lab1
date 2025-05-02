#include<iostream>
#include <algorithm>

#include <string>
#include <stdio.h>


#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"




double square(double x) {return x * x;}

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);

}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
 
class Ray {
    public:
    
        Ray(const Vector& origin, const Vector& direction) : origin(origin), direction(direction) {}
        Vector origin;
        Vector direction;
        
    };

class BoundingBox {
    public:
        Vector Bmin;
        Vector Bmax;
        
        BoundingBox() {
            Bmin = Vector(999999, 999999, 999999);
            Bmax = Vector(-999999, -999999, -999999);
        }
        
        void find_min(const Vector& v) {
            for (int i = 0; i < 3; i++) {
                Bmin[i] = std::min(Bmin[i], v[i]);
                Bmax[i] = std::max(Bmax[i], v[i]);
            }
        }
        
        bool intersect(const Ray& ray, double& t0, double& t1) const {
            t0 = -999999;
            t1 = 999999;
            for (int i = 0; i < 3; ++i) {
                double raydir = 1.0 / ray.direction[i];
                double T0 = (Bmin[i] - ray.origin[i]) * raydir;
                double T1  = (Bmax[i] - ray.origin[i]) * raydir;
        
                if (T0 > T1) std::swap(T0, T1);
        
                t0 = std::max(t0, T0);
                t1 = std::min(t1, T1);
        
                if (t0 > t1) return false; 
            }
            return true;
        }
    };

class Geometry {
    public:
        Vector albedo;
        bool mirror;
        bool transparent;
        double refraction_index;
        
        Geometry(const Vector& albedo , bool mirror = false, bool transparent = false, double refraction_index = 1.5)
            : albedo(albedo), mirror(mirror), transparent(transparent), refraction_index(refraction_index) {}
        
        virtual bool intersect(const Ray& r, double& t, Vector& N, Vector& P) const = 0;
        virtual ~Geometry() {}
    };
        

class Sphere : public Geometry {
public:
    Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false,
        double refraction_index = 1.5, bool reverse_normal = false) :
        Geometry(albedo, mirror, transparent, refraction_index), center(center), radius(radius), reverse_normal(reverse_normal){}

    bool intersect(const Ray& r, double& t, Vector& N, Vector& P) const override{
        Vector oc = r.origin - center;
        double u_dot_oc = dot(r.direction, oc);
        double delta = square(u_dot_oc) - oc.norm2() + square(radius);
    
        if (delta < 0) {
            return false;
        }

        double sqrt_delta = sqrt(delta);
        double u_dot_co = dot(r.direction, center - r.origin);
        
        double t1 = u_dot_co - sqrt_delta;
        double t2 = u_dot_co + sqrt_delta;
        if (t1 > 0) {
            t = t1;
        }
        else if (t2 > 0){
        t = t2;
        }
        else {
        return false;
        }
    
        P = r.origin + t * r.direction;
        N = P - center;
        if (reverse_normal) N = Vector(0,0,0) - N;
        N.normalize();
        return true;
    }
    
    
    


    Vector center;
    double radius;
    bool reverse_normal;
};


class Tree_BVH {
    public:
        BoundingBox bbox;
        int start, end;
        Tree_BVH* left = nullptr;
        Tree_BVH* right = nullptr;
    };
    

   
class TriangleIndices {
    public:
        TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
        };
        int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
        int uvi, uvj, uvk;  // indices within the uv coordinates array
        int ni, nj, nk;  // indices within the normals array
        int group;       // face group
    };
     
     
class TriangleMesh : public Geometry {
    public:
        ~TriangleMesh() {}
        TriangleMesh()
        : Geometry(Vector(1, 1, 1), false, false, 1.5) {}
        Tree_BVH* root = nullptr;

        Vector barycenter(const TriangleIndices& tri) const {
            return (vertices[tri.vtxi] + vertices[tri.vtxj] + vertices[tri.vtxk]) / 3.0;
        }

        Tree_BVH* build_bvh(int start, int end) {
            Tree_BVH* node = new Tree_BVH();
            node->start = start;
            node->end = end;
        
            for (int i = start; i < end; ++i) {
                node->bbox.find_min(vertices[indices[i].vtxi]);
                node->bbox.find_min(vertices[indices[i].vtxj]);
                node->bbox.find_min(vertices[indices[i].vtxk]);
            }
        
        
            Vector diag = node->bbox.Bmax - node->bbox.Bmin;

            int axis;
            if (diag[0] > diag[1] && diag[0] > diag[2]){
                axis = 0;
            } else if (diag[1] > diag[2]){
                axis = 1;
            } else {
                axis = 2;
            }

            double axis_mid = (node->bbox.Bmin[axis] + node->bbox.Bmax[axis]) / 2.0;
        
            int pivot = start;
            for (int i = start; i < end; ++i) {
                if (barycenter(indices[i])[axis] < axis_mid) {
                    std::swap(indices[i], indices[pivot]);
                    pivot++;
                }
            }
        
            if (pivot == start || pivot == end || end - start < 5) return node; 
        
            node->left = build_bvh(start, pivot);
            node->right = build_bvh(pivot, end);
            return node;
        }
        
    

        void readOBJ(const char* obj) {
     
            char matfile[255];
            char grp[255];
     
            FILE* f;
            f = fopen(obj, "r");
            int curGroup = -1;
            while (!feof(f)) {
                char line[255];
                if (!fgets(line, 255, f)) break;
     
                std::string linetrim(line);
                linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
                strcpy(line, linetrim.c_str());
     
                if (line[0] == 'u' && line[1] == 's') {
                    sscanf(line, "usemtl %[^\n]\n", grp);
                    curGroup++;
                }
     
                if (line[0] == 'v' && line[1] == ' ') {
                    Vector vec;
     
                    Vector col;
                    if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                        col[0] = std::min(1., std::max(0., col[0]));
                        col[1] = std::min(1., std::max(0., col[1]));
                        col[2] = std::min(1., std::max(0., col[2]));
     
                        vertices.push_back(vec);
                        vertexcolors.push_back(col);
     
                    } else {
                        sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                        vertices.push_back(vec);
                    }
                }
                if (line[0] == 'v' && line[1] == 'n') {
                    Vector vec;
                    sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    normals.push_back(vec);
                }
                if (line[0] == 'v' && line[1] == 't') {
                    Vector vec;
                    sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                    uvs.push_back(vec);
                }
                if (line[0] == 'f') {
                    TriangleIndices t;
                    int i0, i1, i2, i3;
                    int j0, j1, j2, j3;
                    int k0, k1, k2, k3;
                    int nn;
                    t.group = curGroup;
     
                    char* consumedline = line + 1;
                    int offset;
     
                    nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                    if (nn == 9) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                        if (nn == 6) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                            if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                            if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                            if (nn == 3) {
                                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                                indices.push_back(t);
                            } else {
                                nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                                indices.push_back(t);
                            }
                        }
                    }
     
                    consumedline = consumedline + offset;
     
                    while (true) {
                        if (consumedline[0] == '\n') break;
                        if (consumedline[0] == '\0') break;
                        nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                        TriangleIndices t2;
                        t2.group = curGroup;
                        if (nn == 3) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                            indices.push_back(t2);
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            k2 = k3;
                        } else {
                            nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                                if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                                if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                j2 = j3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                                if (nn == 2) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                    if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                    if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    k2 = k3;
                                    indices.push_back(t2);
                                } else {
                                    nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                    if (nn == 1) {
                                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                        consumedline = consumedline + offset;
                                        i2 = i3;
                                        indices.push_back(t2);
                                    } else {
                                        consumedline = consumedline + 1;
                                    }
                                }
                            }
                        }
                    }
     
                }
     
            }
            fclose(f);
     
        }
     
        std::vector<TriangleIndices> indices;
        std::vector<Vector> vertices;
        std::vector<Vector> normals;
        std::vector<Vector> uvs;
        std::vector<Vector> vertexcolors;
        BoundingBox bbox;

        TriangleMesh(const Vector& albedo,  bool mirror = false, bool transparent = false, double refraction_index = 1.5)
        : Geometry(albedo, mirror, transparent, refraction_index) {}

        
        bool intersect_bvh(const Ray& ray, double& t, Vector& N, Vector& P, Tree_BVH* node) const {
            double t0, t1;
            if (!node->bbox.intersect(ray, t0, t1)) return false;
        
            if (!node->left && !node->right) {
                bool hit = false;
                double t_min = 999999;
        
                for (int i = node->start; i < node->end; ++i) {
                    const Vector& A = vertices[indices[i].vtxi];
                    const Vector& B = vertices[indices[i].vtxj];
                    const Vector& C = vertices[indices[i].vtxk];
        
                    Vector e1 = B - A;
                    Vector e2 = C - A;
                    Vector N = cross(e1, e2); 
                    
                    double denominator = dot(ray.direction, N);
                    if (fabs(denominator) < 1e-6) continue; 
                    
                    
                    
                    double beta = dot(e2, cross(A - ray.origin, ray.direction)) / denominator;
                    double gamma = - dot(e1, cross(A - ray.origin, ray.direction)) / denominator;
                    double alpha = 1.0 - beta - gamma;
                    
                    if (beta < 0 || gamma < 0 || alpha < 0) continue;

            
                    double t_temp =  dot(A - ray.origin, N) / denominator;
                    if (t_temp < 1e-4 || t_temp > t_min) continue;

                    t_min = t_temp;
                    P = ray.origin + t_temp * ray.direction;
                    N = cross(e1, e2);
                    N.normalize();
                    hit = true;
                    
                }
        
                if (hit) {
                    t = t_min;
                    return true;
                }
                return false;
            }
        
            bool hit_left = false, hit_right = false;
            double t_left = 999999, t_right = 999999;
            Vector P_left, N_left, P_right, N_right;
        
            if (node->left) hit_left = intersect_bvh(ray, t_left, N_left, P_left, node->left);
            if (node->right) hit_right = intersect_bvh(ray, t_right, N_right, P_right, node->right);
        
            if (hit_left && (!hit_right || t_left < t_right)) {
                t = t_left; 
                P = P_left;
                N = N_left;
                return true;
            } else if (hit_right) {
                t = t_right;
                P = P_right;
                N = N_right;
                return true;
            }
        
            return false;
        }
        
        bool intersect(const Ray& r, double& t, Vector& N, Vector& P) const override {
            return intersect_bvh(r, t, N, P, root);
        }

        
        
    };



class Scene {
    public:
        
        std::vector<Geometry*> objects;
        Scene(const std::vector<Geometry*>& objects) : objects(objects) {}
        Vector getColor(const Ray& ray, int depth) {
            if (depth < 0) {
            return Vector(0, 0, 0);
            }
            
            double t_min = 9999999;
            Geometry* hit_obj = nullptr;
            Vector P, N;
            
            double t;
            for (auto* obj : objects) {
                Vector N1, P1;
                if (obj->intersect(ray, t, N1, P1) && t < t_min) {
                    t_min = t;
                    N = N1;
                    P = P1;
                    hit_obj = obj;
                }
            }
            
            if (!hit_obj) return Vector(0, 0, 0);

            
            if (!hit_obj->mirror && !hit_obj->transparent) {
                Vector Lo(0,0,0);
                Vector light_pos(-10, 20, 30);
                double I =2e10;
            
                Vector L_minus_P = light_pos - P;
                double d2 = L_minus_P.norm2();
                double d = sqrt(d2);
                Vector light_direction = L_minus_P / d;
            
                Ray shadow_ray(P + 1e-4 * N, light_direction);
                bool light = true; 
            
                for (auto* obj : objects) {
                    double t2;
                    Vector N2, P2;
                    if (obj->intersect(shadow_ray, t2, N2, P2) && t2 < d) {
                        light = false;
                        break;
                    }
                }
            
                if (light) {
                    Vector direct_light = I / (4 * M_PI * d2) * (hit_obj->albedo / M_PI) * std::max(0.0, dot(N, light_direction));
                    return direct_light;
                } else {
                
                Ray random_ray(P + 1e-4 * N, random_cos(N));
                Lo = Lo + hit_obj->albedo * getColor(random_ray, depth - 1);
                return Lo; 
                }
    
            } else {
            
            double n1 = 1.0, n2 = hit_obj->refraction_index;

            double wi_N = dot(ray.direction, N);
            bool exit = wi_N > 0;
            if (exit) {
                std::swap(n1, n2);
                N = Vector(0, 0, 0) - N;
                wi_N = dot(ray.direction, N);
            }
            
            double R0 = square((n1 - n2) / (n1 + n2));
            double R = R0 + (1 - R0) * pow(1 - abs(wi_N), 5);
            double u = (double)rand() / RAND_MAX;
            
            Vector next_dir;
            Vector offset;
            double k = 1 - square(n1/n2) * (1 - wi_N * wi_N);
    
            if (hit_obj->mirror || (hit_obj->transparent && u < R) || k < 0) {
                next_dir = ray.direction - 2 * dot(ray.direction, N) * N;
                offset = 1e-4 * N;
            } else {
                next_dir = n1/n2 * ray.direction - (n1/n2 * wi_N + sqrt(k)) * N;
                offset = -1e-4 * N;
            }
            
            next_dir.normalize();
            Ray next_ray(P + offset, next_dir);
            return hit_obj->albedo * getColor(next_ray, depth - 1); 
            }
        }
    
        Vector random_cos(const Vector &N){
            double r1 = (double)rand() / RAND_MAX;
            double r2 = (double)rand() / RAND_MAX;
            double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
            double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
            double z = sqrt(r2);
            
            Vector T1;
            if (abs(N[0]) <= abs(N[1]) && abs(N[0]) <= abs(N[2]))
                T1 = Vector(0, -N[2], N[1]);
            else if (abs(N[1]) <= abs(N[0]) && abs(N[1]) <= abs(N[2]))
                T1 = Vector(-N[2], 0, N[0]);
            else
                T1 = Vector(-N[1], N[0], 0);
            
            T1.normalize();
            Vector T2 = cross(N, T1);
            
            return x * T1 + y * T2 + z * N;
        }
            
            
};

    

    
      