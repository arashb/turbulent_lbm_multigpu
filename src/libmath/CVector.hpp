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


#ifndef __CVECTOR_HH
#define __CVECTOR_HH
/*
 * CHANGELOG:
 *
 * 2008-10-22: martin schreiber
 * 	created this class to handle vectors in CVector<number of elements><type> way
 *
 * 	this should make it easyer to create more flexible filters like gradient
 * 	independend on the vectorsize
 */


/// template nD vector handler
/*
 * it's not possible to create templates dependend on dimension AND type parameters while specifying
 * only the dimension. therefore this fake template has to be created while the methods are implemented
 * directly in the CVector?.hpp header files
 */
template<int N, typename T>
class CVector
{
};

#include "CVector2.hpp"
#include "CVector3.hpp"
#include "CVector4.hpp"

#endif
