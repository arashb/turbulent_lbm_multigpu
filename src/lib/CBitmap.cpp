/*
 * Copyright 2010 Martin Schreiber
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "CBitmap.hpp"
#include <stdio.h>

#include <string>

/**
 * allocate memory for bitmap with given parameters
 */
CBitmap24::CBitmap24(int p_width, int p_height)	:
	width(p_width),
	height(p_height)
{
	data = new char[width*height*3];
}

/**
 * http://en.wikipedia.org/wiki/BMP_file_format
 */
bool CBitmap24::save(std::string &filename)
{
	FILE *file = fopen(filename.c_str(), "w");
	if (!file)
	{
		error << "cannot open file '" << filename.c_str() << "'" << std::endl;
		return false;
	}

	struct
	{
		unsigned short magic_number;
	} file_header_1;

	struct
	{
		unsigned int file_size;
		unsigned short res1;
		unsigned short res2;
		unsigned int offset_data;
	} file_header_2;

	int padding = (4-(3*width)%4)%4;

	file_header_1.magic_number = 'B'+'M'*256;
	file_header_2.file_size = 54+(3*width+padding)*height;
	file_header_2.res1 = 0;
	file_header_2.res2 = 0;
	file_header_2.offset_data = 54;

	size_t x;
	x = fwrite(&file_header_1, sizeof(file_header_1), 1, file);
	x = fwrite(&file_header_2, sizeof(file_header_2), 1, file);

	struct
	{
		unsigned int header_size;
		unsigned int width;
		unsigned int height;
		unsigned short color_planes;
		unsigned short bpp;
		unsigned int compression_method;
		unsigned int image_size;
		unsigned int h_res;
		unsigned int v_res;
		unsigned int colors_in_palette;
		unsigned int important_colors;
	} image_header;

	image_header.header_size = 40;
	image_header.width = width;
	image_header.height = height;
	image_header.color_planes = 1;
	image_header.bpp = 24;
	image_header.compression_method = 0;
	image_header.image_size = width*height*3;
	image_header.h_res = 11811;
	image_header.v_res = 11811;
	image_header.colors_in_palette = 0;
	image_header.important_colors = 0;

	x = fwrite(&image_header, sizeof(image_header), 1, file);

	// The rows in BMPs always contain a multiple of four pixels, so pad at the end with black pixels.
	// The resolution, however, may be arbitrary. That's what the iPad is good for.
	// write data
	if (padding == 0)
	{
		x = fwrite(data, 3*width*height, 1, file);
	}
	else
	{
		for (int y = 0; y < height; y++)
		{
			x = fwrite(&data[y*width*3], 3*width, 1, file);
			x = fwrite(data, padding, 1, file);
		}
	}

	fclose(file);

	return true;
}

/**
 * free bitmap memory
 */
CBitmap24::~CBitmap24()
{
	if (data)
		delete data;
}
