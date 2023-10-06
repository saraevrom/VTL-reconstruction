from __future__ import annotations
from typing import Union, Optional

import numpy as np
import numpy

HALFPI = np.pi/2


def assert_real(x):
    assert isinstance(x, int) or isinstance(x, float)


def assert_real(x):
    assert isinstance(x, int) or isinstance(x, float)




class Vector2(object):
    '''
    Essential class that represents 2D vector
    '''

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def unpack(self):
        return self.x, self.y

    @staticmethod
    def make_vector(item:Union[list[Union[int, float]], Vector2]):
        if isinstance(item, Vector2):
            return item
        if isinstance(item, list):
            return Vector2(item[0],item[1])

    @staticmethod
    def from_list(l):
        if isinstance(l, list):
            return [Vector2(x[0], x[1]) for x in l]
        else:
            return Vector2(l[0], l[1])

    @property
    def xy(self):
        '''
        Swizzling support
        :return:
        '''
        return self

    @property
    def yx(self):
        '''
        Swizzling support
        :return:
        '''
        return Vector2(self.y, self.x)

    @staticmethod
    def zero():
        '''
        same as Vector2(0,0)
        :return:
        '''
        return Vector2(0,0)

    @staticmethod
    def right():
        '''
        same as Vector2(1,0)
        :return:
        '''
        return Vector2(1, 0)

    @staticmethod
    def left():
        '''
        same as Vector2(-1,0)
        :return:
        '''
        return Vector2(-1, 0)

    @staticmethod
    def up():
        '''
        same as Vector2(0,1)
        :return:
        '''
        return Vector2(0, 1)

    @staticmethod
    def down():
        '''
        same as Vector2(0, -1)
        :return:
        '''
        return Vector2(0, -1)

    def __repr__(self):
        return f"V2[{self.x}, {self.y}]"

    def __add__(self, other):
        return Vector2(
            self.x + other.x,
            self.y + other.y,
        )

    def __abs__(self):
        return Vector2(abs(self.x), abs(self.y))

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __sub__(self, other):
        return Vector2(
            self.x - other.x,
            self.y - other.y,
        )

    def __mul__(self, other):
        return Vector2(self.x * other, self.y * other)

    def __rmul__(self, other):
        return Vector2(self.x * other, self.y * other)

    def __truediv__(self, other):
        return Vector2(self.x / other, self.y / other)

    def dot(self, other):
        '''
        Dot product of two vectors
        :param other: second vector
        :return:
        '''
        return self.x * other.x + self.y * other.y

    def length(self):
        '''
        Length of vector
        :return:
        '''
        return (self.dot(self)) ** 0.5

    def normalized(self):
        '''
        Get unit vector with same direction.
        :return:
        '''
        return self / self.length()

    def get_glsl_declaration(self):
        return f"vec2({self.x}, {self.y})"

    def rotate(self, angle):
        '''
        Create vector rotated by angle CCW
        :param angle: angle of rotation
        :return:
        '''
        return Vector2(
            np.cos(angle) * self.x - np.sin(angle) * self.y,
            np.sin(angle) * self.x + np.cos(angle) * self.y
        )

    def __neg__(self):
        return Vector2(-self.x, -self.y)


    def to_numpy(self):
        return np.array([self.x, self.y])

    def to_list(self):
        return [self.x, self.y]

class Vector3(object):
    '''
    Essential class that represents 3D vector
    '''
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def unpack(self):
        return self.x, self.y, self.z

    @staticmethod
    def make_vector(item:Union[list[Union[int, float]], Vector3]):
        if isinstance(item, Vector3):
            return item
        if isinstance(item, list):
            return Vector3(item[0],item[1],item[2])

    @staticmethod
    def from_list(l:list):
        if isinstance(l[0],list):
            return [Vector3(x[0],x[1],x[2]) for x in l]
        else:
            return Vector3(l[0],l[1],l[2])

    @staticmethod
    def zero():
        '''
        Same as Vector3(0, 0, 0)
        :return: Constant vector
        '''
        return Vector3(0, 0, 0)

    @staticmethod
    def right():
        '''
        Same as Vector3(1, 0, 0)
        :return: Constant vector
        '''
        return Vector3(1,0,0)

    @staticmethod
    def left():
        '''
        Same as Vector3(-1, 0, 0)
        :return: Constant vector
        '''
        return Vector3(-1, 0, 0)

    @staticmethod
    def up():
        '''
        Same as Vector3(0, 1, 0)
        :return: Constant vector
        '''
        return Vector3(0, 1, 0)

    @staticmethod
    def down():
        '''
        Same as Vector3(0, -1, 0)
        :return: Constant vector
        '''
        return Vector3(0, -1, 0)

    @staticmethod
    def forward():
        '''
        Same as Vector3(0, 0, -1).
        (Take a note that coordinate system is left-handed)
        :return: Constant vector
        '''
        return Vector3(0, 0, -1)

    @staticmethod
    def back():
        '''
        Same as Vector3(0, 0, 1).
        (Take a note that coordinate system is left-handed)
        :return: Constant vector
        '''
        return Vector3(0, 0, 1)

    def __repr__(self):
        return f"V3[{self.x},{self.y},{self.z}]"

    def __add__(self, other):
        if isinstance(other, Vector3):
            return Vector3(
                self.x + other.x,
                self.y + other.y,
                self.z + other.z,
            )
        else:
            return other.__radd__(self)

    def __abs__(self):
        return Vector3(abs(self.x), abs(self.y), abs(self.z))

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __sub__(self, other):
        return Vector3(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            raise RuntimeError("Cannot multiply vector and quaternion (Did you mean opposite order?)")
        return Vector3(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other):
        return Vector3(self.x * other, self.y * other, self.z * other)

    def __truediv__(self, other):
        return Vector3(self.x / other, self.y / other, self.z / other)

    def dot(self, other):
        '''
        Dot product of two vectors
        :param other: second vector
        :return: dot product
        '''
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        '''
        Cross product of two vectors
        :param other: second vector
        :return: cross product
        '''
        return Vector3(
            x=self.y * other.z - self.z * other.y,
            y=self.z * other.x - self.x * other.z,
            z=self.x * other.y - self.y * other.x,
        )

    def length(self):
        '''
        Get length of vector
        :return:
        '''
        return (self.dot(self)) ** 0.5

    def sqr_len(self):
        return self.dot(self)

    def normalized(self):
        '''
        Get unit vector with same direction.
        :return:
        '''
        l = self.length()
        return self / l

    def get_glsl_declaration(self):
        return f"vec3({self.x}, {self.y}, {self.z})"

    def __neg__(self):
        return Vector3(-self.x, -self.y, -self.z)

    def to_numpy(self):
        return np.array([self.x, self.y, self.z])

    def to_list(self):
        return [self.x, self.y, self.z]


# QUATERNION MULTIPLICATION RULES
# ROW * COLUMN
#    1  i  j  k
# 1  1  i  j  k
# i  i -1  k -j
# j  j -k -1  i
# k  k  j -i -1


class Quaternion(object):
    '''
    Essential class that represents rotation.
    '''
    def __init__(self, w, x, y, z):
        # assert_real(w)
        # assert_real(x)
        # assert_real(y)
        # assert_real(z)
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Q[{self.w}, V({self.x}, {self.y}, {self.z})]"

    @staticmethod
    def identity():
        '''
        Identity quaternion. Has no rotation effect
        :return:
        '''
        return Quaternion(1, 0, 0, 0)

    @staticmethod
    def from_axis_rotation(angle: float, axis: Vector3, backend=numpy):
        '''
        Quaternion that applies rotation around axis with given angle
        :param angle: angle of rotation
        :param axis: axis of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.from_scalar_vector(backend.cos(angle / 2), axis.normalized()*backend.sin(angle / 2))

    @staticmethod
    def from_to_rotation(dir1: Vector3, dir2: Vector3, hint_dir:Optional[Vector3]=None):
        '''
        Quaternion that will rotate vector `dir1` to `dir2`
        :param dir1: initial vector
        :param dir2: final vector
        :return:
        '''

        dir1 = dir1.normalized()
        dir2 = dir2.normalized()

        full_cos = dir1.dot(dir2)
        if (2 * (1 + full_cos)) == 0:
            if hint_dir is None:
                # Vector got flipped. Let's guess axis.
                if abs(dir2.dot(Vector3.back())) == 1:
                    # vector is aligned with z axis
                    axis = Vector3.right()
                else:
                    x = -dir2.y
                    y = dir2.x
                    axis = Vector3(x, y, 0).normalized()
            else:
                axis = hint_dir.normalized()
            # Pure vector quaternion
            return Quaternion.from_vector(axis)
        else:
            axis = dir1.cross(dir2)
            full_sin = axis.length()
            half_cos = ((1 + full_cos) / 2) ** 0.5
            half_sin = full_sin / (2 * (1 + full_cos)) ** 0.5

            axis_norm = axis.normalized()

            return Quaternion.from_scalar_vector(half_cos, axis_norm * half_sin)

    @staticmethod
    def look_at(position: Vector3, target: Vector3, up=Vector3(0, 1, 0), positiveZ=False):
        '''
        Quaternion that makes rotation for object with given positon that makes it look towards target.
        Vector `up` aligns object up.
        :param position: object position
        :param target: target to look at
        :param up: desired up direction for object
        :return:
        '''
        if positiveZ:
            origin = Vector3(0, 0, 1)
        else:
            origin = Vector3(0, 0, -1)
        main_up = Vector3(0.0, 1.0, 0.0)
        up = up.normalized()
        direction = (target - position).normalized()
        up = (up - direction * direction.dot(up)).normalized()
        main_rot = Quaternion.from_to_rotation(origin, direction)
        rotated_main_up = main_rot.rotate_vector(main_up)
        second_rot = Quaternion.from_to_rotation(rotated_main_up, up, hint_dir=direction)
        return second_rot * main_rot

    @staticmethod
    def rotate_xy(angle=HALFPI, backend=numpy):
        '''
        Shortcut for rotation around axis Z.
        :param angle: angle of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.from_axis_rotation(angle, Vector3(0, 0, 1), backend)

    @staticmethod
    def rotate_yx(angle=HALFPI, backend=numpy):
        '''
        Shortcut for rotation around axis -Z.
        :param angle: angle of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.rotate_xy(-angle, backend)

    @staticmethod
    def rotate_yz(angle=HALFPI, backend=numpy):
        '''
        Shortcut for rotation around axis X.
        :param angle: angle of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.from_axis_rotation(angle, Vector3(1, 0, 0), backend)

    @staticmethod
    def rotate_zy(angle=HALFPI, backend=numpy):
        '''
        Shortcut for rotation around axis -X.
        :param angle: angle of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.rotate_yz(-angle,backend)

    @staticmethod
    def rotate_zx(angle=HALFPI, backend=numpy):
        '''
        Shortcut for rotation around axis Y.
        :param angle: angle of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.from_axis_rotation(angle, Vector3(0, 1, 0),backend)

    @staticmethod
    def rotate_xz(angle=HALFPI, backend=numpy):
        '''
        Shortcut for rotation around axis -Y.
        :param angle: angle of rotation
        :param backend: trigonometry backend
        :return:
        '''
        return Quaternion.rotate_zx(-angle, backend)

    @staticmethod
    def from_scalar_vector(w, v:Vector3):
        '''
        Create quaternion from scalar/vector pair.
        :param w: scalar part
        :param v: vector part
        :return:
        '''
        return Quaternion(w, v.x, v.y, v.z)

    @staticmethod
    def from_scalar(w):
        '''
        Create pure scalar quaternion with zero vector part.
        :param w: scalar part
        :return:
        '''
        return Quaternion(w, 0, 0, 0)

    @staticmethod
    def from_vector(v: Vector3):
        '''
        Create pure vector quaternion with zero scalar part.
        :param v: vector part
        :return:
        '''
        return Quaternion(0, v.x, v.y, v.z)

    @staticmethod
    def _to_quaternion(v):
        if isinstance(v, Quaternion):
            return v
        else:
            return Quaternion.from_scalar(v)

    @staticmethod
    def euler(alpha,beta,gamma):
        '''
        Quaternion from euler angles
        :param alpha:
        :param beta:
        :param gamma:
        :return:
        '''
        q1 = Quaternion.rotate_xy(alpha)
        beta_axis = q1.rotate_vector(Vector3.right())
        q2 = Quaternion.from_axis_rotation(beta, beta_axis)
        q12 = q2*q1
        gamma_axis = q2.rotate_vector(Vector3.back())
        q3 = Quaternion.from_axis_rotation(gamma, gamma_axis)
        return q3*q12

    def scalar(self):
        '''
        Get scalar part of quaternion
        :return: scalar part
        '''
        return self.w

    def vector(self):
        '''
        Get vector part of quaternion
        :return: vector part
        '''
        return Vector3(self.x, self.y, self.z)

    def __add__(self, other):
        other = self._to_quaternion(other)
        return Quaternion(
            self.w + other.w,
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __radd__(self, other):
        return Quaternion.from_scalar(other) + self

    def __rsub__(self, other):
        return Quaternion.from_scalar(other) - self

    def __sub__(self, other):
        other = self._to_quaternion(other)
        return Quaternion(
            self.w - other.w,
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def conj(self):
        '''
        Get conjugated quaternion. It performs opposite rotation.
        :return: Quaternion(w, -x, -y, -z)
        '''
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    def dot(self, other):
        '''
        Dot product of two quaternions.
        :param other: second quaternion
        :return: dot product
        '''
        return self.w * other.w + self.x * other.x + self.y * other.y + self.z * other.z

    def __mul__(self, other):
        if isinstance(other,Quaternion):
            scalar1 = self.scalar()
            scalar2 = other.scalar()
            vector1 = self.vector()
            vector2 = other.vector()

            scalar = scalar1 * scalar2 - vector1.dot(vector2)
            #print(scalar1)
            #print(scalar2)
            #print(scalar1*vector2)
            vector = vector2*scalar1 + vector1*scalar2 + vector1.cross(vector2)
            assert isinstance(vector,Vector3)
            return Quaternion.from_scalar_vector(scalar, vector)
        elif isinstance(other,Vector3):
            return self.rotate_vector(other)
        elif isinstance(other, float) or isinstance(other, int):
            return Quaternion(self.w*other, self.x*other, self.y*other, self.z*other)
        else:
            return other.__rmul__(self)

    def __truediv__(self, other):
        return Quaternion(self.w / other, self.x / other, self.y / other, self.z / other)

    def __neg__(self):
        return Quaternion(-self.w, -self.x, -self.y, -self.z)

    def length(self):
        '''
        Get length of quaternion
        :return: length of quaternion
        '''
        return (self.dot(self)) ** 0.5

    def normalized(self):
        '''
        Get quaternion with unit length and same direction.
        :return: unit quaternion
        '''
        return self / self.length()

    def inverse(self):
        '''
        Get inverted quaternion.
        :return: inverted quaternion
        '''
        return self.conj() / self.dot(self)

    def rotate_vector(self, vec: Vector3):
        '''
        Rotates vector.
        :param vec: vector to be rotated.
        :return: rotated vector.
        '''
        rotated = self * Quaternion.from_vector(vec) * self.inverse()
        return rotated.vector()

    def to_numpy_xyzw(self):
        return np.array([self.x, self.y, self.z, self.w])

    def to_numpy_wxyz(self):
        return np.array([self.w, self.x, self.y, self.z])

    def to_list_xyzw(self):
        return [self.x, self.y, self.z, self.w]

    def to_list_wxyz(self):
        return [self.w, self.x, self.y, self.z]
