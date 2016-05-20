#include "incs.h"
#include "lib/EasyBMP.h"

#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

class image_utils {

public:
    ndvector read_in_image(const char* filename,int &height,int &width);
    int number_of_colors() const;
    void write_to_file(const char* filename, ndvector &image,const int &height,const int &width);
    vector colormap(int size);

private:
    BMP bmp_image;
    ndvector vec_image;
    int num_cols;
};

#endif
