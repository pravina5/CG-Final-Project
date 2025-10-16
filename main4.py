from rit_window import *
from cgI_engine import *
from shapes import *
import numpy as np
from PIL import Image


class TransformStack:
    def __init__(self):
        self.stack = []

    def push_matrix(self, matrix):
        if self.stack:
            self.stack.append(np.dot(self.stack[-1], matrix))
        else:
            self.stack.append(matrix)

    def pop_matrix(self):
        if self.stack:
            return self.stack.pop()
        else:
            return None


def draw_cube(im, transform_stack, translation, rotatez, rotatey, rotatex, scale, viewT, projectionT, r, g, b):
    transform_stack.push_matrix(translation @ rotatez @ rotatey @ rotatex @ scale)
    myEngine.drawTrianglesTextures(cube, cube_idx, cube_uv, im, transform_stack.stack[-1], viewT, projectionT)
    transform_stack.pop_matrix()


# def draw_cone(im, transform_stack, translation, rotatez, rotatey, rotatex, scale, viewT, projectionT, r, g, b):
#     transform_stack.push_matrix(translation @ rotatez @ rotatey @ rotatex @ scale)
#     # myEngine.drawTrianglesTextures(cone, cone_idx,cone_uv,im, transform_stack.stack[-1], viewT, projectionT)
#     myEngine.drawTriangles3Dwireframe(cone, cone_idx, transform_stack.stack[-1], viewT, projectionT, r, g, b, 0.0, 0.0,
#                                       0.0)
#
#     transform_stack.pop_matrix()
#

def draw_cone(im, transform_stack, translation, rotatez, rotatey, rotatex, scale, viewT, projectionT, r, g, b):
    transform_stack.push_matrix(translation @ rotatez @ rotatey @ rotatex @ scale)
    myEngine.drawTrianglesCheckerboard(cone, cone_idx, cone_uv, [0.6, 0.4, 0.2], [0.96, 0.96, 0.86], 0.1,
                                       transform_stack.stack[-1], viewT, projectionT)

    transform_stack.pop_matrix()

# Beige or Cream: [0.96, 0.96, 0.86] for a light beige or [1, 0.98, 0.94] for a creamy white.
# Teal or Turquoise: [0, 0.5, 0.5] for a darker teal or [0.25, 0.88, 0.82] for a brighter turquoise.
# Mustard Yellow: [1, 0.86, 0.35] for a rich mustard yellow.

def draw_cylinder(im, transform_stack, translation, rotatez, rotatey, rotatex, scale, viewT, projectionT, r, g, b):
    transform_stack.push_matrix(translation @ rotatez @ rotatey @ rotatex @ scale)
    myEngine.drawTrianglesTextures(cylinder, cylinder_idx, cylinder_uv, im, transform_stack.stack[-1], viewT,
                                   projectionT)
    transform_stack.pop_matrix()


def default_action():
    # clear the FB
    myEngine.clearFB(0, 0, 0.0)
    myEngine.defineViewWindow(700, 0, 700, 0)

    # define your viewing and projection transforms (frustum) here
    eye = glm.vec3(0.0, 2.0, -4.0)
    lookat = glm.vec3(0.0, 0.0, -10.0)
    up = glm.vec3(0.0, 1.0, 0.0)

    viewT = myEngine.lookAt(eye, lookat, up)
    projectionT = myEngine.frustum3D(-6, 6, -6, 6, 8, 15)
    im = Image.open("tiledstones.jpg")

    # draw your Scene here
    transform_stack = TransformStack()

    # cube
    transform_stack.push_matrix(myEngine.translate3D(0.5, 1.0, -10.0))
    draw_cube(im, transform_stack, myEngine.translate3D(0, -1.3, 0), myEngine.rotateZ(0.0), myEngine.rotateY(20.0),
              myEngine.rotateX(10.0), myEngine.scale3D(3.8, 3.24, 2.8), viewT, projectionT, 1, 0, 0)
    transform_stack.pop_matrix()

    # cylinder left back
    transform_stack.push_matrix(myEngine.translate3D(-2.1, -0.1, -10.0))

    draw_cylinder(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
                  myEngine.rotateX(0.0), myEngine.scale3D(1.0, 3.0, 0.5), viewT, projectionT, 1, 1, 0)

    transform_stack.pop_matrix()

    # cylinder left front
    transform_stack.push_matrix(myEngine.translate3D(-1.0, -0.1, -7.5))

    draw_cylinder(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
                  myEngine.rotateX(0.0), myEngine.scale3D(1.0, 3.0, 0.5), viewT, projectionT, 1, 1, 0)

    transform_stack.pop_matrix()

    # cylinder right front

    transform_stack.push_matrix(myEngine.translate3D(2.0, -0.1, -8.0))

    draw_cylinder(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
                  myEngine.rotateX(0.0), myEngine.scale3D(1.0, 3.0, 0.5), viewT, projectionT, 1, 1, 0)

    transform_stack.pop_matrix()

    # cylinder right back

    transform_stack.push_matrix(myEngine.translate3D(2.6, -0.1, -12.0))

    draw_cylinder(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
                  myEngine.rotateX(0.0), myEngine.scale3D(1.0, 3.0, 0.5), viewT, projectionT, 1, 1, 0)

    transform_stack.pop_matrix()

    # cone front right
    transform_stack.push_matrix(myEngine.translate3D(1.9, 1.9, -7.8))

    draw_cone(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
              myEngine.rotateX(0.0), myEngine.scale3D(1.0, 1.0, 0), viewT, projectionT, 0.6, 0.4, 0.2)

    transform_stack.pop_matrix()

    # cone back right
    transform_stack.push_matrix(myEngine.translate3D(2.6, 1.9, -12.0))

    draw_cone(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
              myEngine.rotateX(0.0), myEngine.scale3D(1.0, 1.0, 0), viewT, projectionT, 0.6, 0.4, 0.2)

    transform_stack.pop_matrix()

    # cone left front

    transform_stack.push_matrix(myEngine.translate3D(-0.9, 1.9, -7.3))

    draw_cone(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
              myEngine.rotateX(0.0), myEngine.scale3D(0.95, 1.0, 0), viewT, projectionT, 0.6, 0.4, 0.2)

    transform_stack.pop_matrix()

    # cone left back

    transform_stack.push_matrix(myEngine.translate3D(-2.1, 1.9, -10.0))

    draw_cone(im, transform_stack, myEngine.translate3D(0, 0.15, 0), myEngine.rotateZ(0.0), myEngine.rotateY(0.0),
              myEngine.rotateX(0.0), myEngine.scale3D(1.0, 1.0, 0.0), viewT, projectionT, 0.6, 0.4, 0.2)

    transform_stack.pop_matrix()

    # Phong shaded moon
    transform_stack.push_matrix(myEngine.translate3D(4.0, 5.0, -16.0))

    myEngine.drawTrianglesPhong(sphere, sphere_idx, sphere_normals, transform_stack.stack[-1], viewT, projectionT,
                                [1.0, 1.0, 1.0], [1.0, 0.8, 0.0], [0.2, 0.4, 0.4], 10.0, [-2.0, 3.0, -2.0],
                                [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], False)


window = RitWindow(800, 800)
myEngine = CGIengine(window, default_action)


def main():
    window.run(myEngine)


if __name__ == "__main__":
    main()
