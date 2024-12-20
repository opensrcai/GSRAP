import numpy as np
from common.camera import Camera
from common.math import normalize


def sample_in_unit_sphere(u: np.ndarray, v: np.ndarray, w: np.ndarray):
    cos_theta = -2 * u + 1
    sin_theta = np.sqrt(1 - cos_theta * cos_theta)
    phi = 2 * np.pi * v
    r = np.power(w, 1.0 / 3.0)

    x = r * sin_theta * np.cos(phi)
    y = r * sin_theta * np.sin(phi)
    z = r * cos_theta

    return np.stack([x, y, z], axis=-1)


def sample_in_unit_sphere_surface(u: np.ndarray, v: np.ndarray):
    cos_theta = -2 * u + 1
    sin_theta = np.sqrt(1 - cos_theta * cos_theta)
    phi = 2 * np.pi * v

    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta

    return np.stack([x, y, z], axis=-1)


def create_spherical_point_cloud(
    center: np.ndarray, radius: float, num_points: int, seed: int, dtype=np.float64
) -> np.ndarray:
    np.random.seed(seed)

    u = np.random.rand(num_points)
    v = np.random.rand(num_points)
    w = np.random.rand(num_points)

    points = sample_in_unit_sphere(u, v, w)
    points = radius * points + center

    return points


def create_cuboid_point_cloud(
    center: np.ndarray, x: float, y: float, z: float, num_points: int, seed: int
):
    np.random.seed(seed)

    # ret = np.zeros((num_point, 3))

    random_values = np.random.rand(num_points, 3) - 0.5

    return center + random_values * np.array([x, y, z])


def to_bearing_vectors(points: np.ndarray, camera: Camera):
    assert points.ndim == 2
    assert points.shape[1] == 3

    points_c = (camera.rotation_cg[None, :, :] @ points[:, :, None]).squeeze(
        -1
    ) + camera.translation_cg[None, :].squeeze(0)

    return normalize(points_c)


def rotate_bearing_vectors(axis: np.ndarray, theta: float, phi: float):
    axis = normalize(axis)


def create_cameras(dtype=np.float64):
    # Set optional axis and up at the global coordinates.
    # Variable optional_axis_g is equal to the z axis of the coordinates of
    # camera.

    camera = Camera(dtype=dtype)

    camera_optical_axis_g = normalize(np.array([-1.0, 1.0, 1.0], dtype=dtype))
    camera_up_g = np.array([0.0, -1.0, 0.0], dtype=dtype)

    # Compute x and y axes of the coordinates of camera.
    camera_xg = normalize(np.cross(camera_up_g, camera_optical_axis_g))
    camera_yg = normalize(np.cross(camera_optical_axis_g, camera_xg))

    # Create global to camera rotation matrix.
    camera.rotation_cg[0, :] = camera_xg
    camera.rotation_cg[1, :] = camera_yg
    camera.rotation_cg[2, :] = camera_optical_axis_g

    # Set camera to global translation.
    # The distance from camera to camera1 is 1.
    camera_translation_gc = np.array(
        [1.0 / np.sqrt(2.0), -1.0 / np.sqrt(2.0), 0.0], dtype=dtype
    )

    # Compute global to camera translation.
    camera.translation_cg = -(camera.rotation_cg.transpose() @ camera_translation_gc)

    return camera
