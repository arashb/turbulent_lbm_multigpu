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


#ifndef FILE_H
#define FILE_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;

class CFile
{
// read out the contents of a file, store the file length to *fileLength and return the pointer to the data
// if terminateNull is true, the content string is terminated by appending simply a 0 byte
//	static void *fileContents(const char *filepath, int *fileLength = 0, bool terminateNull = false);

public:

//static void *fileContents(const char *filepath, int *fileLength, bool terminateNull)
/*
 * read the contents to string
 * \param filePath	filepath to source file
 * \param data		data where to store the filecontents
 * \param errorLog	error output
 */
static bool fileContents(const string &filePath, string &p_data, string &errorLog)	//bool terminateNull)
{
	int length;
	char *data;

	FILE *file = fopen(filePath.c_str(), "rb");
	if (file == NULL)
	{
		perror("fileContent() fopen");
		return false;
	}

	fseek(file, 0, SEEK_END);
	length = ftell(file);
	fseek(file, 0, SEEK_SET);

	data = new char[length];

	if (!data)
	{
		errorLog ="fileContent() out of memory!";
		return false;
	}

	int readLength;
	readLength = fread(data, 1, length, file);
	fclose(file);

	if (readLength != length)
	{
		errorLog = "readLength != fileLength";
		return false;
	}

	p_data.assign(data, readLength);
	return true;
}

static bool storeToFile(const char *filepath, void *data, size_t length)
{
	FILE *file = fopen(filepath, "w");
	if (file == NULL)
	{
		perror("fileContent() fopen");
		return 0;
	}

	// TODO: error check
	size_t x = fwrite(data, length, 1, file);
	fclose(file);

	if (x < length)	return false;

	return true;
}


/**
 * load data from file
 * \param filepath	file to load
 * \param data		storage area
 * \param max_length	maximum length
 * \return		-1: error	else: loaded bytes
 */
static int loadFromFile(const char *filepath, void *data, int max_length)
{
	int length;

	FILE *file = fopen(filepath, "w");
	if (file == NULL)
	{
		perror("fileContent() fopen");
		return 0;
	}

	fseek(file, 0, SEEK_END);
	length = ftell(file);
	fseek(file, 0, SEEK_SET);

	if (max_length != -1)
		if (length > max_length-1)
		{
			cerr << "not enough space (filesize: " << length << endl;
			return -1;
		}

	// TODO: error check
	int read_length = fread(data, length, 1, file);
	fclose(file);

	return read_length;
}


};

#endif
