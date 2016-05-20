#include <mpi.h>
#include <iostream>
#include <cmath>
#include "image_utils.h"
#include "incs.h"

const int MAX_COLOR = 256;
int height,width;

void denoise_image(image_utils &util,ndvector &image,int iter,double eta,double beta,int rank, int nprocs);
val find_color_value(ndvector &colormap,ndvector &image,ndvector &hidden_nodes,int i,int j,double eta,double beta,bool BW,int number_of_colors);
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

	MPI_Init(&argc,&argv);

    image_utils util;
    ndvector image;
	int rank,nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	if(rank == 0) {
		//only read in image once
		image = util.read_in_image(input_image,height,width);
	}

	MPI_Bcast(&height,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&width,1,MPI_INT,0,MPI_COMM_WORLD);

	if(rank != 0)
		image = ndvector(height*width);

	//send image vector to all processes
	MPI_Bcast(&image.front(),height*width,MPI_INT,0,MPI_COMM_WORLD);	

	if(rank == 0) {
		std::cerr << "Preparing to denoise..." << std::endl;
	}
	denoise_image(util,image,iterations,eta,beta,rank,nprocs);

	if(rank == 0)
		std::cerr << "Finished denoising, writing to file..." << std::endl;
	
	//only write to file once
	if(rank == 0)
		util.write_to_file(output_image,image,height,width);

	MPI_Finalize();
	

    return 0;
}

void denoise_image(image_utils &util,ndvector &image,int iterations,double eta,double beta,int rank,int nprocs) {
	
    bool BW = false;
    vector colormap;
	int number_of_colors;

	if(rank == 0) {
		colormap = util.colormap(height*width);
		number_of_colors = util.number_of_colors();
	}

	MPI_Bcast(&number_of_colors,1,MPI_INT,0,MPI_COMM_WORLD);

	if(rank != 0)
		colormap = vector(number_of_colors);
	
	MPI_Bcast(&colormap.front(),number_of_colors,MPI_INT,0,MPI_COMM_WORLD);

    if(number_of_colors == 2) {
		if(rank == 0)
			BW = true;
		MPI_Bcast(&BW,1,MPI::BOOL,0,MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
		std::cerr << "Reading processes..." << std::endl;
	vector recvcounts = vector(nprocs);
	vector displs = vector(nprocs);
	vector start = vector(nprocs);
	vector end = vector(nprocs);
	// recv = # of rows per process
	int recv = height/nprocs;
	for(int i=0; i < nprocs; ++i) {
		recvcounts[i] = recv;
	}
	for(int i=0; i < recv*nprocs - height; ++i) {
		//distribute remainder over processes
		++recvcounts[i];
	}
	int total=0;
	int starttot = 0;
	for(int i=0; i < nprocs; ++i) {
		start[i] = starttot;
		end[i] = starttot = starttot + recvcounts[i];
		
		//*width to get total # of pixels
		recvcounts[i] *= width;
		
		displs[i] = total;
		total += recvcounts[i];
	}

	if(rank == 0)
		std::cerr << "Starting denoising iterations..." << std::endl;

    ndvector hidden_nodes = ndvector(image.size());
	
    #pragma omp parallel for
	for(int i=0;i<height*width;++i)
		hidden_nodes[i] = image[i];

	for(int iter=0; iter<iterations; ++iter) {
        ndvector new_hidden = ndvector(image.size());
		//each process only iterates over its rows
#pragma omp parallel for
        for(int i=start[rank]; i < end[rank]; ++i) {
			//each row is thread-parallelized 
            //#pragma omp parallel for num_threads(8)
            for(int j=0; j < width; ++j) {
                val color = find_color_value(colormap,image,hidden_nodes,i,j,eta,beta,BW,number_of_colors);
                new_hidden[i*width+j] = color;
            }
        }

#pragma omp parallel for
		for(int i=start[rank];i<end[rank]; ++i) {
			for(int j=0;j<width;++j) {
				hidden_nodes[i*width+j] = new_hidden[i*width+j];
			}
		}
		//hidden_nodes = ndvector(new_hidden);
		//combine all hidden_nodes vectors so everyone has all info for next iteration
		MPI_Allgatherv(MPI_IN_PLACE,0,MPI_INT,&hidden_nodes.front(),recvcounts.data(),displs.data(),MPI_INT,MPI_COMM_WORLD);
    }

	#pragma omp parallel for
	for(int i=0; i < height*width; ++i)
		image[i] = hidden_nodes[i];

}

val find_color_value(vector &colormap,ndvector &image,ndvector &hidden_nodes,int i,int j,double eta,double beta,bool BW,int number_of_colors) {

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
