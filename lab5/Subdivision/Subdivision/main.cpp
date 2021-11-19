#include <string.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include "halfedge.h"
#include "catmullclark.h"

int option;

std::vector<Mesh*> meshes;
int g_iteration = 0;
std::string g_filename;

std::vector<Vector3f> draw_vertices;
std::vector<Vector3f> draw_normals;

int g_vsize, g_nsize, g_tsize, g_isize;

void catmull_clark();
void catmull_clark_subdivision();

// only used for the initial mesh, normals for subdivision are calculated there
void calculate_normals() {
    draw_normals.clear();
    int inc = meshes[g_iteration]->quad ? 4 : 3;
    for (int i = 0; i < draw_vertices.size(); i+=inc) {
        // each three will make a face
        Vector3f a = draw_vertices[i];
        Vector3f b = draw_vertices[i+1];
        Vector3f c = draw_vertices[i+2];
        Vector3f d;

        if (meshes[g_iteration]->quad) {
            d = draw_vertices[i+3];

        }
            // find for each vertex
        Vector3f ab = a - b;
        Vector3f bc = b - c;
        Vector3f normal;
        cross(&normal, &ab, &bc);
        normal.normalize();
        draw_normals.push_back(normal);
        draw_normals.push_back(normal);
        draw_normals.push_back(normal);

        if (meshes[g_iteration]->quad) {
            draw_normals.push_back(normal);
        }
    }
}

void print_draw() {
    int inc = meshes[g_iteration]->quad ? 4 : 3;
    for (int i = 0; i < draw_vertices.size(); i+=inc) {
        printf("v ");
        draw_vertices[i].print();
        printf("v ");
        draw_vertices[i+1].print();
        printf("v ");
        draw_vertices[i+2].print();

        if (inc == 4) {
            draw_vertices[i+3].print();
            //draw_normals[i+3].print();
        }
        printf("\n");
    }
}

void print_indices() {
    int inc = meshes[g_iteration]->quad ? 4 : 3;
    for (int i = 0; i < draw_vertices.size(); i+=inc) {
        printf("f %d %d %d\n", i+1, i+2, i+3);
    }
}

void prepare_to_draw(int it) {
    draw_vertices.clear();
    draw_normals.clear();
    if (it < 0 || it > g_iteration) {return;}
    for (std::vector<Face*>::iterator f = meshes[it]->glfaces.begin();
         f != meshes[it]->glfaces.end(); f++) {
        HalfEdge *e = (*f)->edge;
        do {
            draw_vertices.push_back(e->start->pos);
            draw_normals.push_back(e->start->normal);
            e = e->next;
        } while (e != (*f)->edge);
    }
}

// Catmull-Clark!!!!!!!!
void catmull_clark_subdivision() {
    generate_face_points(meshes[g_iteration], meshes[g_iteration-1]); // 1
    generate_edge_points(meshes[g_iteration], meshes[g_iteration-1]); // 2
    generate_new_vertices(meshes[g_iteration], meshes[g_iteration-1]); // 3
    connect_new_mesh(meshes[g_iteration], meshes[g_iteration-1]); // 4
}

void catmull_clark() {
    Mesh *mesh = new Mesh(); // create a new mesh
    mesh->quad = true; // catmull generates quads
    meshes.push_back(mesh); // add the new mesh
    g_iteration++;

    printf("-------------------start %d--------------------\n", g_iteration);
    catmull_clark_subdivision(); // perform the subdivision
	printf("-------------------done %d-------------------\n", g_iteration);

    prepare_to_draw(g_iteration); // copy mesh vertices into draw_arrays
}

void init_data() {
    g_iteration = 0;
    Mesh *mesh = new Mesh();
    meshes.push_back(mesh);

    if (!meshes[g_iteration]->load_file(g_filename)) { // load data
        printf("Mesh file couldn't be laoded\n");
        exit(0);
    }
    // create the half edge data structure and load the draw array
    init_mesh(meshes[g_iteration], draw_vertices);

    calculate_normals();
}

int main(int argc, char** argv) {
    /*
     if (!(argc == 2)) {
     printf("usage: ./hw5 option \n");
     menu();
     return 0;
     } */
    int option = 2;
    if (option == 0) {
        g_filename = "input/quad.obj";
    }
    else if (option == 1) {
        g_filename = "input/triangle.obj";
    }
    else if (option == 2) {
        g_filename = "input/venus.obj";
    }
    else if (option == 3) { // cow
        g_filename = "input/cow.obj";
    }
    else if (option == 4) { // bunny
        g_filename = "input/bunny.obj";
    }
    else if (option == 5) { // teddy
        g_filename = "input/teddy.obj";
    }
    else if (option == 6) {
        g_filename = "input/stell.obj";
    }
    else if (option == 7) {
        g_filename = "input/ico.obj";
    }
    else if (option == 8) {
        g_filename = "input/torus.obj";
    }

    // initial setup
    init_data();
	catmull_clark();
	system("pause");

    return 0;
}
