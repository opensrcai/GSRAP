import math
import multiprocessing
from typing import cast

import gsrap
import numpy as np
from common.utils import (
    Camera,
    create_cameras,
    create_spherical_point_cloud,
    to_bearing_vectors,
)


def main():
    dtype = np.float64
    num_points = 500

    center = np.array(
        [0.5 / np.sqrt(2.0), -0.5 / np.sqrt(2.0), 0.5 / np.sqrt(2.0)], dtype=dtype
    )  # Center of spherical point cloud.

    radius: float = 0.5 / math.sqrt(2.0)

    points = create_spherical_point_cloud(
        center, radius, num_points, 123456789, dtype=dtype
    )

    camera1 = Camera(dtype)
    camera2 = create_cameras(dtype)

    bearings1 = to_bearing_vectors(points, camera1)
    bearings2 = to_bearing_vectors(points, camera2)

    matches = np.arange(num_points, dtype=np.uint32)[:, None].repeat(2, axis=1)

    essential_solver_policy = gsrap.EssentialSolverPolicy()

    essential_solver_policy.flags = (
        gsrap.EssentialSolverPolicyFlags.ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE
    )
    essential_solver_policy.ransac_policy.flags = (
        gsrap.RansacPolicyFlags.RANSAC_POLICY_FLAGS_EARLY_STOP
        | gsrap.RansacPolicyFlags.RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE
    )

    essential_solver_policy.ransac_policy.probability = 0.999999
    essential_solver_policy.ransac_policy.num_threads = multiprocessing.cpu_count()
    essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 30
    essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 100

    essential_solver_policy.inlier_thr = 1e-12

    ret = gsrap.compute_essential_matrix(
        essential_solver_policy, bearings1, bearings2, matches
    )

    essential_solver_result = ret[0]
    if essential_solver_result is not None:
        essential_solver_result = cast(
            gsrap.EssentialSolverResult, essential_solver_result
        )

        print("R:\n", essential_solver_result.rotation)
        print("det R: ", np.linalg.det(essential_solver_result.rotation))
        print("gt R:\n", camera2.rotation_cg)
        print("---")
        print("t:\n", essential_solver_result.translation)
        print("gt t:\n", camera2.translation_cg)


if __name__ == "__main__":
    main()
