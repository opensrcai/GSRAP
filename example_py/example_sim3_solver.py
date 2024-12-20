import multiprocessing
from typing import cast

import gsrap
import numpy as np
from common.utils import create_cuboid_point_cloud
from scipy.spatial.transform import Rotation


def main():
    num_points = 5000

    center = np.array([0, 0, 0])  # Center of spherical point cloud.
    edge_length = 1.0  # Edge length of cuboid point cloud

    gt_scale = 1.25

    axis = np.array([1, 1, 1], dtype=np.double)
    axis /= np.linalg.norm(axis, ord=2)
    angle = 0.25 * np.pi

    gt_r = Rotation.from_rotvec(angle * axis).as_matrix()

    gt_t = np.array([1.0, 0.25, 0.5])

    points1 = create_cuboid_point_cloud(
        center, edge_length, edge_length, edge_length, num_points, 123456
    )

    points2 = (gt_scale * gt_r @ points1.transpose()).transpose() + gt_t

    matches = np.arange(num_points, dtype=np.uint32)[:, None].repeat(2, axis=1)

    sim3_solver_policy = gsrap.Sim3SolverPolicy()

    sim3_solver_policy.flags = int(
        gsrap.Sim3SolverPolicyFlags.SIM3_SOLVER_POLICY_FLAGS_NONE
    ) | int(gsrap.Sim3SolverPolicyFlags.SIM3_SOLVER_POLICY_FLAGS_REFINE)

    sim3_solver_policy.ransac_policy.flags = int(
        gsrap.RansacPolicyFlags.RANSAC_POLICY_FLAGS_EARLY_STOP
    ) | int(
        gsrap.RansacPolicyFlags.RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE
    )
    sim3_solver_policy.ransac_policy.num_threads = multiprocessing.cpu_count()
    sim3_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10
    sim3_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000

    sim3_solver_policy.inlier_thr = 1e-12
    gsrap.compute_sim3_transformation(sim3_solver_policy, points1, points2, matches)
    ret = gsrap.compute_sim3_transformation(
        sim3_solver_policy, points1, points2, matches
    )

    sim3_solver_result = ret[0]
    if sim3_solver_result is not None:
        sim3_solver_result = cast(gsrap.Sim3SolverResult, sim3_solver_result)
        print("scale:\n", sim3_solver_result.scale, "\n")
        print("gt scale:\n", gt_scale, "\n")

        print("R :", sim3_solver_result.rotation, "\n")
        print("gt R :", gt_r, "\n")

        print("t :", sim3_solver_result.translation, "\n")
        print("gt t :", gt_t, "\n")

    print("RANSAC interation: ", ret[1].num_iteration)


if __name__ == "__main__":
    main()
