def normalize(x):
    import numpy as np
    norm = np.linalg.norm(x)
    if np.isclose(norm, 1.0):
        return x
    if np.isclose(norm, 0.0):
        raise ValueError(f"Zero-length vector detected for the vector.")
    return x / norm



def robust_align(a, b):
    """
    Portable function: Perfectly align two vectors. Handle 180 degree case.
    Return rotation matrix
        - a: rotate vector
        - b: reference vector
    """
    import numpy as np
    from scipy.spatial.transform import Rotation
    # Normalize two vectors
    norm_a, norm_b = np.linalg.norm(a), np.linalg.norm(b)
    if (norm_a < 1e-8) or (norm_b < 1e-8):
        print(norm_a, norm_b)
        raise ValueError(f"Zero-length vector detected for one of the vector.")
    a = a / norm_a
    b = b / norm_b

    dot = np.dot(a, b)

    # 0deg case: consider as aligned
    if np.isclose(dot, 1.0):
        return np.eye(3) # identity matrix

    # 180° case: flip
    if np.isclose(dot, -1.0):
        axis = np.array([1,0,0]) # x-axis
        if abs(a[0]) > 0.9: # if vector is too close to rotation axis, the rotation operation is unstable
            axis = np.array([0,1,0]) # y-axis
        rot_axis = np.cross(a, axis) # To get an axis perp to a -> rotation axis
        rot_axis /= np.linalg.norm(rot_axis)
        return Rotation.from_rotvec(np.pi * rot_axis).as_matrix() # rotate 180° = inverse = - indentity matrix

    rot_axis = np.cross(a, b)
    angle = np.arccos(np.clip(dot, -1, 1))

    # cos(x) = cos(-x), but in cross product sin(x) = -sin(-x) => Combine two, return robust align
    rot = Rotation.from_rotvec(angle * rot_axis / np.linalg.norm(rot_axis)) # rotate angle° matrix
    return rot.as_matrix() # rotation matrix/operator


