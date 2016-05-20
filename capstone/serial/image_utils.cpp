#include "image_utils.h"

ndvector image_utils::read_in_image(const char* filename,int &height,int &width) {

    bmp_image.ReadFromFile(filename);
    height = bmp_image.TellHeight(), width = bmp_image.TellWidth();
    vec_image = ndvector(height*width);

    for(int i=0; i < height; ++i) {
        for(int j=0; j < width; ++j) {
            vec_image[i*width+j] = bmp_image(j,i)->Red;
        }
    }


    return vec_image;
}

int image_utils::number_of_colors() const {
    return num_cols;
}

ndvector image_utils::colormap(int size) {
    std::set<val> colors;
    for(int i=0; i < size; ++i)
		colors.insert(vec_image[i]);
	num_cols = colors.size();
    ndvector color_vec = ndvector(colors.size());
	int i=0;
	for(auto &c : colors) {
		color_vec[i] = c;
		++i;
	}
    return color_vec;
}

void image_utils::write_to_file(const char* filename,ndvector &image,const int &height, const int &width) {
    for(int i=0; i < height; ++i) {
        for(int j=0; j < width; ++j) {
            bmp_image(j,i)->Red = image[i*width+j];
            bmp_image(j,i)->Blue = image[i*width+j];
            bmp_image(j,i)->Green = image[i*width+j];
        }
    }
    bmp_image.WriteToFile(filename);
}


