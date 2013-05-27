/*
 * common.h
 *
 *  Created on: May 16, 2013
 *      Author: Arash Bakhtiari
 */

#ifndef COMMON_H_
#define COMMON_H_

#define FLAG_OBSTACLE	(1 << 0)
#define FLAG_FLUID	(1 << 1)
#define FLAG_VELOCITY_INJECTION	(1 << 2)
#define FLAG_GHOST_LAYER (1 << 3)
#define FLAG_GHOST_LAYER_BETA (FLAG_GHOST_LAYER | (1 << 4))

extern CVector<3,int> E0 	;
extern CVector<3,int> E1	;
extern CVector<3,int> E2	;
extern CVector<3,int> E3	;

extern CVector<3,int> E4	;
extern CVector<3,int> E5	;
extern CVector<3,int> E6	;
extern CVector<3,int> E7	;

extern CVector<3,int> E8	;
extern CVector<3,int> E9	;
extern CVector<3,int> E10	;
extern CVector<3,int> E11	;

extern CVector<3,int> E12	;
extern CVector<3,int> E13	;
extern CVector<3,int> E14	;
extern CVector<3,int> E15	;

extern CVector<3,int> E16	;
extern CVector<3,int> E17	;
extern CVector<3,int> E18	;

extern CVector<3,int> lbm_units[];


#endif /* COMMON_H_ */