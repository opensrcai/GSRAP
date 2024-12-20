import math
import multiprocessing
from typing import cast

import gsrap
import numpy as np
from common.utils import (
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

    camera = create_cameras()

    bearings = to_bearing_vectors(points, camera)

    matches = np.arange(num_points, dtype=np.uint32)[:, None].repeat(2, axis=1)

    # print("center: ", center)
    # print("radius: ", radius)
    # print("camera:\n", camera.rotation_cg, "\n---\n", camera.translation_cg)
    # print("points: ", points.shape, " ", points.dtype)
    # print("bearings: ", bearings.shape, " ", bearings.dtype)
    # print("matches: ", matches.shape, " ", matches.dtype)

    pnp_solver_policy = gsrap.PnpSolverPolicy()
    # pnp_solver_policy.ransac_policy   = {};

    pnp_solver_policy.flags = (
        gsrap.PnpSolverPolicyFlags.PNP_SOLVER_POLICY_FLAGS_NONE
        | gsrap.PnpSolverPolicyFlags.PNP_SOLVER_POLICY_FLAGS_REFINE
    )

    pnp_solver_policy.ransac_policy.flags = (
        gsrap.RansacPolicyFlags.RANSAC_POLICY_FLAGS_EARLY_STOP
        | gsrap.RansacPolicyFlags.RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE
    )

    pnp_solver_policy.ransac_policy.num_threads = multiprocessing.cpu_count()
    pnp_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10
    pnp_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000

    pnp_inlier_check_params = gsrap.PnpInlierCheckParamsUsingBearingVector()
    pnp_inlier_check_params.inlier_thr = 1e-6
    pnp_solver_policy.pnp_inlier_check_params = pnp_inlier_check_params

    ret = gsrap.solve_pnp_problem(pnp_solver_policy, bearings, points, matches)

    pnp_solver_result = ret[0]
    if pnp_solver_result is not None:
        pnp_solver_result = cast(gsrap.PnpSolverResult, pnp_solver_result)

        print("R:\n", pnp_solver_result.rotation)
        print("det R: ", np.linalg.det(pnp_solver_result.rotation))
        print("gt R:\n", camera.rotation_cg)
        print("---")
        print("t:\n", pnp_solver_result.translation)
        print("gt t:\n", camera.translation_cg)


if __name__ == "__main__":
    main()
