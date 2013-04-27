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


#ifndef CBITMAP_HPP
#define CBITMAP_HPP

#include <string>
#include "CError.hpp"


/**
 * \brief	store image data to a bitmap file
 *
 * store the 24 bit rgb image data in the array data[] to a 24 bit bitmap file
 */
class CBitmap24
{
public:
	CError error;		///< error handler

	int width;			///< bitmap width
	int height;			///< bitmap height
	char *data;			///< allocated bitmap data where the 24 bit bitmap data has to be loaded before the bitmap can be stored

	CBitmap24();
	CBitmap24(int p_width, int p_height);
	~CBitmap24();

	bool save(std::string &file);
};

#endif
