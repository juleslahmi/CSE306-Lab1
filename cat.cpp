#include<iostream>
#include <algorithm>
#include <omp.h>
#include <chrono>


#include "classes.cpp"

 
int main() {
    auto start = std::chrono::high_resolution_clock::now();
    int W = 512;
    int H = 512;
    double angle = 60*M_PI/180;
    double gamma = 2.2;
    
    int number_of_rays = 1000;
    

    Vector camera_origin(0, 0, 60);

    std::vector<Geometry*> objects;

    TriangleMesh* mesh = new TriangleMesh();
    mesh->readOBJ("cadnav.com_model/Models_F0202A090/cat.obj");
    for (Vector& v : mesh->vertices) {
        v = 0.6 * v + Vector(0, -10, 0);
    }


    mesh->root = mesh->build_bvh(0, mesh->indices.size());

    objects.push_back(new Sphere(Vector(0, -1000, 0), 990, Vector(0.3, 0.4, 0.7)));
    objects.push_back(new Sphere(Vector(0, 1000, 0), 940, Vector(0.2, 0.5, 0.9)));
    objects.push_back(new Sphere(Vector(0, 0, 1000), 940, Vector(0.9, 0.4, 0.3)));
    objects.push_back(new Sphere(Vector(0, 0, -1000), 940, Vector(0.4, 0.8, 0.7)));
    objects.push_back(new Sphere(Vector(1000, 0, 0), 950, Vector(0.9, 0.2, 0.9)));
    objects.push_back(new Sphere(Vector(-1000, 0, 0), 950, Vector(0.6, 0.5, 0.1)));
    objects.push_back(mesh);
    
    Scene scene(objects);


    Vector light(-10, 20, 40);
    
    double z = -W/(2 * tan(angle/2));

    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1) 
    for (int i = 0; i < H; i++) {
        std::cout << i << std::endl;
        for (int j = 0; j < W; j++) {
            
            Vector pixelColor( 0. , 0. , 0. );
            double x = j - W/2 + 0.5; 
            double y = H/2 - i - 0.5;
            Vector n(x,y,z);
            n.normalize();

            Ray r(camera_origin, n );
            Vector color_sum(0, 0, 0);

            for (int s = 0; s < number_of_rays; ++s) {
                double u1 = (double)rand() / RAND_MAX;
                double u2 = (double)rand() / RAND_MAX;
                double dx = 0.5 * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2); 
                double dy = 0.5 * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    
                double x = j - W / 2.0 + 0.5 + dx ;
                double y = H / 2.0 - i - 0.5 + dy;
    
                Vector dir(x, y, z);
                dir.normalize();
                Ray ray(camera_origin, dir);
    
                color_sum = color_sum + scene.getColor(ray, 8);
            }
            Vector color = color_sum / number_of_rays;
                
            Vector color_gamma(
                pow(color[0], 1.0/gamma),
                pow(color[1], 1.0/gamma),
                pow(color[2], 1.0/gamma)
            );
                
            //std::cout << color_gamma.data[0] << " " << color_gamma.data[1] << " " << color_gamma.data[2] << std::endl;
            if (color_gamma[0] > 255) {
                color_gamma[0] = 255;
            }
            if (color_gamma[1] > 255) {
                color_gamma[1] = 255;
            }
            if (color_gamma[2] > 255) {
                color_gamma[2] = 255;
            }
            image[(i * W + j) * 3 + 0] = color_gamma.data[0];
            image[(i * W + j) * 3 + 1] = color_gamma.data[1];
            image[(i * W + j) * 3 + 2] = color_gamma.data[2];

        }
        }
    stbi_write_png("cat2.png", W, H, 3, &image[0], 0);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}


