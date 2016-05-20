#include <mpi.h>
#include <iostream>
#include <cmath>
#include "image_utils.h"
#include "incs.h"

const int MAX_COLOR = 256;
int height,width;

void denoise_image(image_utils &util,ndvector &image,int iter,double eta,double beta);
val find_color_value(ndvector &colormap,ndvector &image,ndvector &hidden_nodes,int i,int j,double eta,double beta,bool BW);
double node_energy(double x,double y,bool BW,double param);

int main(int argc,char **argv) {
    if(argc < 4) {
        std::cout << "Usage: ./image_denoise [input image] [output image] [num iterations] [eta] [beta]" << std::endl;
        return 1;
    }
    const char *input_image = argv[1];
    const char *output_image = argv[2];
    int iterations = atoi(argv[3]);
    double eta = 2.1, beta = 1.0;
    if(argc > 4)
        eta = atof(argv[4]);
    if(argc > 5)
        beta = atof(argv[5]);

    image_utils util;
    ndvector image;
	image = util.read_in_image(input_image,height,width);

	std::cerr << "Preparing to denoise..." << std::endl;
	denoise_image(util,image,iterations,eta,beta,rank,nprocs);

	std::cerr << "Finished denoising, writing to file..." << std::endl;
	
	util.write_to_file(output_image,image,height,width);

    return 0;
}

void denoise_image(image_utils &util,ndvector &image,int iterations,double eta,double beta) {
	
    bool BW = false;
    vector colormap;
	int number_of_colors;

	colormap = util.colormap(height*width);
	number_of_colors = util.number_of_colors();

    if(number_of_colors == 2) {
		BW = true;
	}

	std::cerr << "Starting denoising iterations..." << std::endl;

    ndvector hidden_nodes = ndvector(image.size());
	
    #pragma omp parallel for
	for(int i=0;i<height*width;++i)
		hidden_nodes[i] = image[i];

	for(int iter=0; iter<iterations; ++iter) {
        ndvector new_hidden = ndvector(image.size());
		//each process only iterates over its rows
#pragma omp parallel for
        for(int i=0; i < height; ++i) {
            for(int j=0; j < width; ++j) {
                val color = find_color_value(colormap,image,hidden_nodes,i,j,eta,beta,BW,number_of_colors);
                new_hidden[i*width+j] = color;
            }
        }


			#pragma omp parallel for
		for(int i=0;i<height; ++i) {
			for(int j=0;j<width;++j) {
				hidden_nodes[i*width+j] = new_hidden[i*width+j];
			}
		}
    }

	#pragma omp parallel for
	for(int i=0; i < height*width; ++i)
		image[i] = hidden_nodes[i];

}

val find_color_value(vector &colormap,ndvector &image,ndvector &hidden_nodes,int i,int j,double eta,double beta,bool BW) {

	double min_energy = 0.0;
	val color_value = -1;
	for(auto &c : colormap) {
		double hidden_observed_energy = 0.0,hidden_hidden_energy = 0.0;
		val new_color_value = c;
		hidden_observed_energy = node_energy(new_color_value,image[i*width+j],BW,eta);
		if(i > 0) {
			hidden_hidden_energy += node_energy(new_color_value,hidden_nodes[(i-1)*width+j],BW,beta);
		}
		if(i < height-1) {
			hidden_hidden_energy += node_energy(new_color_value,hidden_nodes[(i+1)*width+j],BW,beta);
		}
		if(j > 0) {
			hidden_hidden_energy += node_energy(new_color_value,hidden_nodes[i*width+j-1],BW,beta);
		}
		if(j < width-1) {
			hidden_hidden_energy += node_energy(new_color_value,hidden_nodes[i*width+j+1],BW,beta);
		}
		double new_energy = hidden_observed_energy + hidden_hidden_energy;
		if(color_value == -1 || new_energy < min_energy) {
			color_value = new_color_value;
			min_energy = new_energy;
		}
	}

	return color_value;
}

double node_energy(double x,double y,bool BW,double param) {
    if(BW) {
        if(x == y)
            return -1.0*param;
        else
            return param;
    } else {
        return param * log(abs(x-y)+1);
    }

}
