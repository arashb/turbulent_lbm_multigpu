#include "src/cl_programs/lbm_defaults.h"

/**
 * COPY KERNEL
 *
 * Copy a block of  data from src to dst.
 * Based on the assumption that src and dst values are saved in x, y, z order.
 * Since each row of block in x direction is copied with one work item,
 * the number of work items to use with this kernel should be (block_size_y*block_size_z)
 *
 */
__kernel void copy_buffer_rect(
		__global T*  	 src,
		__const			 int src_offset,
		__const 		 int src_origin_x,
		__const  		 int src_origin_y,
		__const          int src_origin_z,
		__const          int src_size_x,
		__const          int src_size_y,
		__const          int src_size_z,
		__global T*  	 dst,
		__const			 int dst_offset,
		__const 		 int dst_origin_x,
		__const          int dst_origin_y,
		__const          int dst_origin_z,
		__const          int dst_size_x,
		__const          int dst_size_y,
		__const          int dst_size_z,
		__const          int block_size_x
		)
{

	int gid_j = get_global_id(0);
	int gid_k = get_global_id(1);

	int src_slice_z_size = src_size_x * src_size_y;
	int dst_slice_z_size = dst_size_x * dst_size_y;

	// index of the origin of the block in the src array
	int src_origin_idx = src_offset + src_origin_x +
						(src_origin_y + gid_j) * src_size_x +
						(src_origin_z + gid_k) * src_slice_z_size;

	// index of the origin of the block in the dst array
	int dst_origin_idx = dst_offset + dst_origin_x +
						(dst_origin_y + gid_j) * dst_size_x +
						(dst_origin_z + gid_k) * dst_slice_z_size;
	int x_idx;
	for	(x_idx = 0; x_idx < block_size_x; x_idx++)
		dst[ dst_origin_idx + x_idx] = src[src_origin_idx + x_idx];

}
