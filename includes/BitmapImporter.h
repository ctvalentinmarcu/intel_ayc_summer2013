/**
 * \file BitmapImporter.h  
 * \brief Contains declarations used to parse and hold the .bmp data.
 * \author Cédric Andreolli (Intel)
 * \date 10 April 2013
 */
#ifndef BITMAP_IMPORTER_H
#define BITMAP_IMPORTER_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <stdint.h>

namespace Accelerate{

/*!
 * \struct HeaderStr
 * \brief This structure holds the .bmp header informations
 * 
 * The HeaderStr structure contains:<br />
 * <ul>
 *	<li>The size of the file</li>
 *	<li>The offset of the data</li>
 *	<li>The width of the image</li>
 *	<li>The height of the image</li>
 *	<li>The number of bits per pixel</li>
 *	<li>Etc.</li>
 * </ul>
 * You don't have to use this structure and you can write/use your own \
 * structure if you want.
 */
typedef struct __attribute__ ((__packed__)){
	//BMP header
	uint16_t magic_number;
	uint32_t size;
	uint32_t reserved;
	uint32_t offset;
	//DIB header
	uint32_t dibSize;
	uint32_t width;
	uint32_t height;
	uint16_t plane;
	uint16_t bit_per_pixel;
	uint32_t compression;
	uint32_t data_size;
	uint32_t hor_res;
	uint32_t vert_res;
	uint32_t color_number;
	uint32_t important;
}HeaderStr;

/*!
 * \struct PixelStr
 * \brief This structure holds a pixel
 * 
 * In the program, a pixel is represented by 3 elements that can define its color.
 * <ul>
 *	<li><b>r: </b>The red value of the pixel</li>
 *	<li><b>g: </b>The green value of the pixel</li>
 *	<li><b>b: </b>The blue value of the pixel</li>
 * </ul>
 */
typedef struct __attribute__ ((__packed__)){
	unsigned char r;
	unsigned char g;
	unsigned char b;
}PixelStr;


/*!
 * \class Image
 * \brief The Image class holds the data of the image. 
 *
 * To instantiate an image from a bitmap, you can use the static function:<br />
 * Image::create_image_from_bitmap(const std::string file_name)<br />
 * 
 */
class Image {
protected:
	PixelStr* pixel_data;/*!< The image data */
	unsigned int width;/*!< The width of the image*/
	unsigned int height;/*!< The height of the image*/
	
	/*!
	 * \brief Default constructor of an image
	 *
	 * This constructor doesn't need to be public.
	 *
	 */
	Image();

	/*!
	 * \brief Copy the column represented by the variable column_number at the right place in a scaled image
	 * \param result The scaled image based on this
	 * \param column_number The column indice to copy in the current image
	 * \param scale The scale factor
	 */
	void copy_column(Accelerate::Image& result, unsigned int column_number, unsigned int scale) const;	
public:
	/*!
	 * \brief The destructor
	 */
	~Image();

	/*!
	 * \brief Returns the pixel at the position [i,j] in the image
	 * 
	 * This function helps you to retrieve the value of a single pixel in the image.
	 * As the image is represented by a single dimension array, to retrieve the pixel
	 * at the position [i,j] you must compute i*image_width + j.
	 *
	 * \param i The rown indice
	 * \param j the column indice
	 * \return The PixelStr at position [i,j]
	 */
	PixelStr get_pixel(unsigned int i, unsigned int j) const;

	/*!
	 * \brief Create a new scaled image from the current image
	 * \param scale The scale factor.
	 * \return A new scaled image. 
	 */
	Image scale_image(unsigned int scale) const;


	void set_pixels(PixelStr* data){ this->pixel_data = data; }
	PixelStr* get_pixels(){ return pixel_data; }
	
	void set_width(const unsigned int width){ this->width = width; }
	unsigned int get_width()const { return this->width; }

	void set_height(const unsigned int height){ this->height = height; }
	unsigned int get_height()const { return this->height; }

	/*!
	 * \brief Creates an image from a bitmap file
	 * \param file_name The path of the bitmap file
	 * \return A new image
	 */
	static Image create_image_from_bitmap(const std::string file_name);
};

}

std::ostream& operator<<(std::ostream &o, const Accelerate::PixelStr& p);
std::ostream& operator<<(std::ostream &o, Accelerate::Image& im);
#endif
