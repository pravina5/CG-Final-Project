from rit_window import *
from vertex import *
import glm as glm
from PIL import Image


class CGIengine:
    def __init__(self, myWindow, defaction):
        self.w_width = myWindow.width
        self.w_height = myWindow.height
        self.win = myWindow
        self.keypressed = 1
        self.default_action = defaction
        self.view_transform = None
        self.depth = []
        self.image = []
        # Initilaize depth and image buffers
        for x in range(0, self.w_width):
            self.depth.append([99999999999] * self.w_height)
            self.image.append([glm.vec3(0.0, 0.0, 0.0)] * self.w_height)

    # go is called on every update of the window display loop
    # have your engine draw stuff in the window.
    def go(self):
        if (self.keypressed == 1):
            # default scene
            self.default_action()

        if (self.keypressed == 2):
            # clear the framebuffer
            self.win.clearFB(0, 0, 0)

        # push the window's framebuffer to the window
        self.win.applyFB()

    def keyboard(self, key):
        if (key == '1'):
            self.keypressed = 1
            self.go()
        if (key == '2'):
            self.keypressed = 2
            self.go()

    def clearFB(self, r, g, b):
        for x in range(self.w_width):
            for y in range(self.w_height):
                self.win.set_pixel(x, y, r, g, b)

    # Calculate cross product between vectors (p1-p0) and ((x,y) - p0)
    def calc_cross_product(self, p0, p1, x, y):
        return ((x - p0.x) * (p1.y - p0.y) - (y - p0.y) * (p1.x - p0.x))

    def identity(self):
        return glm.mat3(1)

    def translate(self, x, y):
        pos = glm.vec2(x, y)
        return glm.translate(pos)

    def scale(self, x, y):
        fac = glm.vec2(x, y)
        return glm.scale(fac)

    def rotate(self, angle):
        # glm calculates cosine and sine using radians and direction of +ve angles is clockwise
        angle = glm.radians((-1) * angle)
        return glm.mat3(glm.cos(angle), (-1) * glm.sin(angle), 0,
                        glm.sin(angle), glm.cos(angle), 0,
                        0, 0, 1)

    def normalize(self, t, b, r, l):
        return glm.mat3(2 / (r - l), 0, 0,
                        0, 2 / (t - b), 0,
                        (((-2) * l) / (r - l)) - 1, (((-2) * b) / (t - b)) - 1, 1)

    def defineViewWindow(self, t, b, r, l):
        self.view_transform = glm.mat3((r - l) / 2, 0, 0,
                                       0, (t - b) / 2, 0,
                                       (r + l) / 2, (t + b) / 2, 1)

    def identity3D(self):
        return glm.mat4(1)

    def translate3D(self, x, y, z):
        pos = glm.vec3(x, y, z)
        return glm.translate(pos)

    def scale3D(self, x, y, z):
        fac = glm.vec3(x, y, z)
        return glm.scale(fac)

    def rotateX(self, angle):
        angle = glm.radians((-1) * angle)
        ans = glm.mat4(1, 0, 0, 0,
                       0, glm.cos(angle), glm.sin(angle), 0,
                       0, (-1) * glm.sin(angle), glm.cos(angle), 0,
                       0, 0, 0, 1)
        return glm.transpose(ans)

    def rotateY(self, angle):
        angle = glm.radians((-1) * angle)
        ans = glm.mat4(glm.cos(angle), 0, (-1) * glm.sin(angle), 0,
                       0, 1, 0, 0,
                       glm.sin(angle), 0, glm.cos(angle), 0,
                       0, 0, 0, 1)
        return glm.transpose(ans)

    def rotateZ(self, angle):
        angle = glm.radians((-1) * angle)
        ans = glm.mat4(glm.cos(angle), glm.sin(angle), 0, 0,
                       (-1) * glm.sin(angle), glm.cos(angle), 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1)
        return glm.transpose(ans)

    def ortho3D(self, l, r, b, t, n, f):
        return glm.orthoRH(l, r, b, t, n, f)

    def is_point_inside(self, p, clip_plane_val, plane_num):
        if plane_num == 0:
            return p.x < clip_plane_val
        elif plane_num == 1:
            return p.y < clip_plane_val
        elif plane_num == 2:
            return p.x > clip_plane_val
        else:
            return p.y > clip_plane_val

    def clip_color(self, p0, p1, p):
        u = glm.sqrt((p[0] - p0.x) ** 2 + (p[1] - p0.y) ** 2) / glm.sqrt((p1.x - p0.x) ** 2 + (p1.y - p0.y) ** 2)
        output = []
        output.append(((1 - u) * p0.r) + (u * p1.r))
        output.append(((1 - u) * p0.g) + (u * p1.g))
        output.append(((1 - u) * p0.b) + (u * p1.b))
        return output

    def compute_intersection(self, p0, p1, clip_plane_val, plane_num):
        if plane_num == 1 or plane_num == 3:
            intersect_coord = [self.clip_y_coord(p0, p1, clip_plane_val), clip_plane_val]
            intersect_color = self.clip_color(p0, p1, intersect_coord)
            return Vertex(intersect_coord[0], intersect_coord[1], 0, None, None, intersect_color[0], intersect_color[1],
                          intersect_color[2])
        else:
            intersect_coord = [clip_plane_val, self.clip_x_coord(p0, p1, clip_plane_val)]
            intersect_color = self.clip_color(p0, p1, intersect_coord)
            return Vertex(intersect_coord[0], intersect_coord[1], 0, None, None, intersect_color[0], intersect_color[1],
                          intersect_color[2])

    def poly_to_triangles(self, vertices):
        output = []
        for i in range(1, len(vertices) - 1):
            output.append(vertices[0])
            output.append(vertices[i])
            output.append(vertices[i + 1])
        return output

    def compute_outcode(self, p, top, bottom, right, left):
        outcode = 0
        # 0001 - left - 1
        # 0010 - right - 2
        # 1000 - up - 8
        # 0100 - down - 4
        if p.x < left:
            outcode |= 1
        if p.x > right:
            outcode |= 2
        if p.y > top:
            outcode |= 8
        if p.y < bottom:
            outcode |= 4
        return outcode

    def clip_x_coord(self, p0, p1, x_val):
        return ((p1.y - p0.y) / (p1.x - p0.x)) * (x_val - p0.x) + p0.y

    def clip_y_coord(self, p0, p1, y_val):
        return ((y_val - p0.y) / (p1.y - p0.y)) * (p1.x - p0.x) + p0.x

    def rasterizeTriangleWireframe(self, p0, p1, p2, wr, wg, wb):
        min_x = min(p0.x, p1.x, p2.x)
        min_y = min(p0.y, p1.y, p2.y)
        max_x = max(p0.x, p1.x, p2.x)
        max_y = max(p0.y, p1.y, p2.y)
        # Area of triangle using cross product
        t = self.calc_cross_product(p0, p1, p2.x, p2.y)
        if t == 0:
            # Early exit in case triangle has zero area (all 3 vettices are colinear)
            return
        # Iterate over extent of triangle
        for x in range(int(min_x), int(max_x + 1)):
            for y in range(int(min_y), int(max_y + 1)):
                t2 = self.calc_cross_product(p0, p1, x, y)
                t0 = self.calc_cross_product(p1, p2, x, y)
                t1 = self.calc_cross_product(p2, p0, x, y)
                # Edge tests to verify if point lies inside the triangle (boundry inclusive)
                is_inside = (t0 >= 0) and (t1 >= 0) and (t2 >= 0)
                # Skip this pixel if it does not lie in the triangle
                if is_inside == False:
                    continue
                # Calculate barycentric coordinates
                c0 = t0 / t
                c1 = t1 / t
                c2 = t2 / t
                # # Interpolate colors
                r = (c2 * p2.r) + (c0 * p0.r) + (c1 * p1.r)
                g = (c2 * p2.g) + (c0 * p0.g) + (c1 * p1.g)
                b = (c2 * p2.b) + (c0 * p0.b) + (c1 * p1.b)
                # Set the pixel color
                z_value = (c2 * p2.z) + (c0 * p0.z) + (c1 * p1.z)
                # Near zero tolerance for rendering edges of wirefrme
                temp = 0.032
                if z_value < self.depth[x][y]:
                    if (c0 <= temp or c1 <= temp or c2 <= temp or c0 == 1.0 or c1 == 1.0 or c2 == 1.0):
                        self.image[x][y] = glm.vec3(wr, wg, wb)
                    else:
                        self.image[x][y] = glm.vec3(r, g, b)
                    self.depth[x][y] = z_value
                self.win.set_pixel(x, y, self.image[x][y].x, self.image[x][y].y, self.image[x][y].z)

    def clipLine(self, P0, P1, top, bottom, right, left):
        while True:
            p0_outcode = self.compute_outcode(P0, top, bottom, right, left)
            p1_outcode = self.compute_outcode(P1, top, bottom, right, left)
            # Trivial Accept
            if p1_outcode | p0_outcode == 0:
                return [P0, P1]
            # Trivial Reject
            if p0_outcode & p1_outcode != 0:
                return []
            # P0 is outside
            if p0_outcode != 0:
                if p0_outcode & 1:
                    # Left
                    P0.y = self.clip_x_coord(P0, P1, left)
                    P0.x = left
                elif p0_outcode & 2:
                    # Right
                    P0.y = self.clip_x_coord(P0, P1, right)
                    P0.x = right
                elif p0_outcode & 4:
                    # Down
                    P0.x = self.clip_y_coord(P0, P1, bottom)
                    P0.y = bottom
                elif p0_outcode & 8:
                    # Up
                    P0.x = self.clip_y_coord(P0, P1, top)
                    P0.y = top
            # P1 is outside
            else:
                if p1_outcode & 1:
                    # Left
                    P1.y = self.clip_x_coord(P1, P0, left)
                    P1.x = left
                elif p1_outcode & 2:
                    # Right
                    P1.y = self.clip_x_coord(P1, P0, right)
                    P1.x = right
                elif p1_outcode & 4:
                    # Down
                    P1.x = self.clip_y_coord(P1, P0, bottom)
                    P1.y = bottom
                elif p1_outcode & 8:
                    # Up
                    P1.x = self.clip_y_coord(P1, P0, top)
                    P1.y = top

    def clipPoly(self, vertices, top, bottom, right, left):
        clip_planes = [right, top, left, bottom]
        output = []
        for index in range(0, len(clip_planes)):
            if len(vertices) == 0:
                continue
            pred = vertices[-1]
            for i in range(0, len(vertices)):
                curr = vertices[i]
                if self.is_point_inside(curr, clip_planes[index], index):
                    if not self.is_point_inside(pred, clip_planes[index], index):
                        clipped = self.compute_intersection(curr, pred, clip_planes[index], index)
                        output.append(clipped)
                    output.append(curr)
                else:
                    if self.is_point_inside(pred, clip_planes[index], index):
                        clipped = self.compute_intersection(pred, curr, clip_planes[index], index)
                        output.append(clipped)
                pred = curr
            vertices = output
            output = []
        return self.poly_to_triangles(vertices)

    def drawTriangles3D(self, vertex_pos, indices, modelT, viewT, projectionT, r, g, b):
        i = 0
        poly_verts = []
        while i < len(indices):
            # Index verts and colors according to specified stride and offsets
            index = indices[i]
            p0 = glm.vec3(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2])
            color0 = [r, g, b]
            i += 1
            index = indices[i]
            p1 = glm.vec3(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2])
            color1 = [r, g, b]
            i += 1
            index = indices[i]
            p2 = glm.vec3(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2])
            color2 = [r, g, b]
            i += 1
            screen_vert_pos0 = projectionT * viewT * modelT * p0
            screen_vert_pos1 = projectionT * viewT * modelT * p1
            screen_vert_pos2 = projectionT * viewT * modelT * p2
            vert0 = Vertex(screen_vert_pos0[0], screen_vert_pos0[1], screen_vert_pos0[2], None, None, color0[0],
                           color0[1],
                           color0[2])
            vert1 = Vertex(screen_vert_pos1[0], screen_vert_pos1[1], screen_vert_pos1[2], None, None, color1[0],
                           color1[1],
                           color1[2])
            vert2 = Vertex(screen_vert_pos2[0], screen_vert_pos2[1], screen_vert_pos2[2], None, None, color2[0],
                           color2[1],
                           color2[2])
            poly = [vert0, vert1, vert2]
            # Clipping
            clipped = self.clipPoly(poly, 1.0, -1.0, 1.0, -1.0)
            poly_verts.extend(clipped)
        nverts = len(poly_verts)
        if nverts < 3:
            return
        if (nverts % 3) != 0:
            print("Bad number of verticies to define a set of triangles ", nverts, " Must be a multiple of 3")
            return
        # simply go through sets of 3 vertices and do raterization
        startidx = 0
        while startidx < nverts:
            p0 = glm.vec3(poly_verts[startidx].x, poly_verts[startidx].y, 1)
            color0 = [poly_verts[startidx].r, poly_verts[startidx].g, poly_verts[startidx].b]
            sp0 = self.view_transform * p0
            poly_verts[startidx] = Vertex(int(sp0[0]), int(sp0[1]), poly_verts[startidx].z, None, None, color0[0],
                                          color0[1],
                                          color0[2])

            p1 = glm.vec3(poly_verts[startidx + 1].x, poly_verts[startidx + 1].y, 1)
            color1 = [poly_verts[startidx + 1].r, poly_verts[startidx + 1].g, poly_verts[startidx + 1].b]
            sp1 = self.view_transform * p1
            poly_verts[startidx + 1] = Vertex(int(sp1[0]), int(sp1[1]), poly_verts[startidx + 1].z, None, None,
                                              color1[0],
                                              color1[1], color1[2])

            p2 = glm.vec3(poly_verts[startidx + 2].x, poly_verts[startidx + 2].y, 1)
            color2 = [poly_verts[startidx + 2].r, poly_verts[startidx + 2].g, poly_verts[startidx + 2].b]
            sp2 = self.view_transform * p2
            poly_verts[startidx + 2] = Vertex(int(sp2[0]), int(sp2[1]), poly_verts[startidx + 2].z, None, None,
                                              color2[0],
                                              color2[1], color2[2])

            self.rasterizeTriangle(poly_verts[startidx], poly_verts[startidx + 1], poly_verts[startidx + 2])
            startidx = startidx + 3

    def drawTriangles3Dwireframe(self, vertex_pos, indices, modelT, viewT, projectionT, r, g, b, wr, wg, wb):
        i = 0
        poly_verts = []
        while i < len(indices):
            # Index verts and colors according to specified stride and offsets
            index = indices[i]
            p0 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            color0 = [r, g, b]
            i += 1
            index = indices[i]
            p1 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            color1 = [r, g, b]
            i += 1
            index = indices[i]
            p2 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            color2 = [r, g, b]
            i += 1
            screen_vert_pos0 = projectionT * viewT * modelT * p0
            screen_vert_pos1 = projectionT * viewT * modelT * p1
            screen_vert_pos2 = projectionT * viewT * modelT * p2
            vert0 = Vertex(screen_vert_pos0[0] / screen_vert_pos0[3], screen_vert_pos0[1] / screen_vert_pos0[3],
                           screen_vert_pos0[2] / screen_vert_pos0[3], None, None, color0[0], color0[1], color0[2])
            vert1 = Vertex(screen_vert_pos1[0] / screen_vert_pos1[3], screen_vert_pos1[1] / screen_vert_pos1[3],
                           screen_vert_pos1[2] / screen_vert_pos1[3], None, None, color1[0], color1[1], color1[2])
            vert2 = Vertex(screen_vert_pos2[0] / screen_vert_pos2[3], screen_vert_pos2[1] / screen_vert_pos2[3],
                           screen_vert_pos2[2] / screen_vert_pos2[3], None, None, color2[0], color2[1], color2[2])
            poly = [vert0, vert1, vert2]
            # Clipping
            clipped = self.clipPoly(poly, 1.0, -1.0, 1.0, -1.0)
            poly_verts.extend(clipped)
        nverts = len(poly_verts)
        if nverts < 3:
            return
        if (nverts % 3) != 0:
            print("Bad number of verticies to define a set of triangles ", nverts, " Must be a multiple of 3")
            return
        # simply go through sets of 3 vertices and do raterization
        startidx = 0
        while startidx < nverts:
            p0 = glm.vec3(poly_verts[startidx].x, poly_verts[startidx].y, 1)
            color0 = [poly_verts[startidx].r, poly_verts[startidx].g, poly_verts[startidx].b]
            sp0 = self.view_transform * p0
            poly_verts[startidx] = Vertex(int(sp0[0]), int(sp0[1]), poly_verts[startidx].z, None, None, color0[0],
                                          color0[1],
                                          color0[2])

            p1 = glm.vec3(poly_verts[startidx + 1].x, poly_verts[startidx + 1].y, 1)
            color1 = [poly_verts[startidx + 1].r, poly_verts[startidx + 1].g, poly_verts[startidx + 1].b]
            sp1 = self.view_transform * p1
            poly_verts[startidx + 1] = Vertex(int(sp1[0]), int(sp1[1]), poly_verts[startidx + 1].z, None, None,
                                              color1[0],
                                              color1[1], color1[2])

            p2 = glm.vec3(poly_verts[startidx + 2].x, poly_verts[startidx + 2].y, 1)
            color2 = [poly_verts[startidx + 2].r, poly_verts[startidx + 2].g, poly_verts[startidx + 2].b]
            sp2 = self.view_transform * p2
            poly_verts[startidx + 2] = Vertex(int(sp2[0]), int(sp2[1]), poly_verts[startidx + 2].z, None, None,
                                              color2[0],
                                              color2[1], color2[2])
            self.rasterizeTriangleWireframe(poly_verts[startidx], poly_verts[startidx + 1], poly_verts[startidx + 2],
                                            wr, wg, wb)
            startidx = startidx + 3

    def lookAt(self, eye, lookat, up):
        return glm.lookAtRH(eye, lookat, up)

    def frustum3D(self, l, r, b, t, n, f):
        return glm.frustumRH_NO(l, r, b, t, n, f)

    def drawTrianglesTextures(self, vertex_pos, indices, uv, im, modelT, viewT, projectionT):
        texture = im  # Load texture image
        tex_width, tex_height = texture.size

        i = 0
        poly_verts = []
        while i < len(indices):
            # Index verts and colors according to specified stride and offsets
            index = indices[i]
            p0 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            u0, v0 = uv[index * 2], uv[index * 2 + 1]
            i += 1
            index = indices[i]
            p1 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            u1, v1 = uv[index * 2], uv[index * 2 + 1]
            i += 1
            index = indices[i]
            p2 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            u2, v2 = uv[index * 2], uv[index * 2 + 1]
            i += 1
            screen_vert_pos0 = projectionT * viewT * modelT * p0
            screen_vert_pos1 = projectionT * viewT * modelT * p1
            screen_vert_pos2 = projectionT * viewT * modelT * p2
            vert0 = Vertex(screen_vert_pos0[0] / screen_vert_pos0[3], screen_vert_pos0[1] / screen_vert_pos0[3],
                           screen_vert_pos0[2] / screen_vert_pos0[3], u0, v0, None, None, None)
            vert1 = Vertex(screen_vert_pos1[0] / screen_vert_pos1[3], screen_vert_pos1[1] / screen_vert_pos1[3],
                           screen_vert_pos1[2] / screen_vert_pos1[3], u1, v1, None, None, None)
            vert2 = Vertex(screen_vert_pos2[0] / screen_vert_pos2[3], screen_vert_pos2[1] / screen_vert_pos2[3],
                           screen_vert_pos2[2] / screen_vert_pos2[3], u2, v2, None, None, None)
            poly = [vert0, vert1, vert2]
            # Clipping
            clipped = self.clipPoly(poly, 1.0, -1.0, 1.0, -1.0)
            poly_verts.extend(clipped)
        nverts = len(poly_verts)
        if nverts < 3:
            return
        if (nverts % 3) != 0:
            print("Bad number of verticies to define a set of triangles ", nverts, " Must be a multiple of 3")
            return
        # simply go through sets of 3 vertices and do raterization
        startidx = 0
        while startidx < nverts:
            p0 = glm.vec3(poly_verts[startidx].x, poly_verts[startidx].y, 1)
            # color0 = [poly_verts[startidx].r, poly_verts[startidx].g, poly_verts[startidx].b]
            uv0 = [poly_verts[startidx].u, poly_verts[startidx].v]
            sp0 = self.view_transform * p0
            poly_verts[startidx] = Vertex(int(sp0[0]), int(sp0[1]), poly_verts[startidx].z, uv0[0], uv0[1], None, None,
                                          None)

            p1 = glm.vec3(poly_verts[startidx + 1].x, poly_verts[startidx + 1].y, 1)
            # color1 = [poly_verts[startidx + 1].r, poly_verts[startidx + 1].g, poly_verts[startidx + 1].b]
            uv1 = [poly_verts[startidx + 1].u, poly_verts[startidx + 1].v]
            sp1 = self.view_transform * p1
            poly_verts[startidx + 1] = Vertex(int(sp1[0]), int(sp1[1]), poly_verts[startidx + 1].z, uv1[0],
                                              uv1[1], None, None, None)

            p2 = glm.vec3(poly_verts[startidx + 2].x, poly_verts[startidx + 2].y, 1)
            # color2 = [poly_verts[startidx + 2].r, poly_verts[startidx + 2].g, poly_verts[startidx + 2].b]
            uv2 = [poly_verts[startidx + 2].u, poly_verts[startidx + 2].v]
            sp2 = self.view_transform * p2
            poly_verts[startidx + 2] = Vertex(int(sp2[0]), int(sp2[1]), poly_verts[startidx + 2].z, uv2[0],
                                              uv2[1], None, None, None)
            self.rasterizeTriangleWithTexture(poly_verts[startidx], poly_verts[startidx + 1], poly_verts[startidx + 2],
                                              texture)
            startidx = startidx + 3

    def rasterizeTriangleWithTexture(self, p0, p1, p2, texture):
        min_x = min(p0.x, p1.x, p2.x)
        min_y = min(p0.y, p1.y, p2.y)
        max_x = max(p0.x, p1.x, p2.x)
        max_y = max(p0.y, p1.y, p2.y)
        # Area of triangle using cross product
        t = self.calc_cross_product(p0, p1, p2.x, p2.y)
        if t == 0:
            # Early exit in case triangle has zero area (all 3 vettices are colinear)
            return
        # Iterate over extent of triangle
        for x in range(int(min_x), int(max_x + 1)):
            for y in range(int(min_y), int(max_y + 1)):
                t2 = self.calc_cross_product(p0, p1, x, y)
                t0 = self.calc_cross_product(p1, p2, x, y)
                t1 = self.calc_cross_product(p2, p0, x, y)
                # Edge tests to verify if point lies inside the triangle (boundry inclusive)
                is_inside = (t0 >= 0) and (t1 >= 0) and (t2 >= 0)
                # Skip this pixel if it does not lie in the triangle
                if is_inside == False:
                    continue
                # Calculate barycentric coordinates
                c0 = t0 / t
                c1 = t1 / t
                c2 = t2 / t
                # Interpolate colors
                u = (c0 * p0.u) + (c1 * p1.u) + (c2 * p2.u)
                v = (c0 * p0.v) + (c1 * p1.v) + (c2 * p2.v)
                print("Pixel at (", x, ",", y, ") - UV coordinates: (", u, ",", v, ")")

                # r = (c2 * p2.r) + (c0 * p0.r) + (c1 * p1.r)
                # g = (c2 * p2.g) + (c0 * p0.g) + (c1 * p1.g)
                # b = (c2 * p2.b) + (c0 * p0.b) + (c1 * p1.b)

                # Sample texture using UV coordinates
                texture_color = texture.getpixel((int(u * (texture.width - 1)), int(v * (texture.height - 1))))

                # Set the pixel color
                z_value = (c2 * p2.z) + (c0 * p0.z) + (c1 * p1.z)
                if z_value < self.depth[x][y]:
                    self.image[x][y] = glm.vec3(texture_color[0] / 255.0, texture_color[1] / 255.0,
                                                texture_color[2] / 255.0)
                    self.depth[x][y] = z_value
                self.win.set_pixel(x, y, self.image[x][y].x, self.image[x][y].y, self.image[x][y].z)

    def calc_lighting(self, vertex_pos, light_pos, normal, k, ocolor, scolor, exponent, light_color,
                      amb_color):
        amb_color = glm.vec3(amb_color[0], amb_color[1], amb_color[2])
        vertex_pos = glm.vec3(vertex_pos[0], vertex_pos[1], vertex_pos[2])
        light_color = glm.vec3(light_color[0], light_color[1], light_color[2])
        ka = k[0]
        kd = k[1]
        ks = k[2]
        view_pos = glm.vec3(0.0, 0.0, 0.0)
        # ambient
        ambient = glm.vec3(ka * amb_color[0], ka * amb_color[1], ka * amb_color[2])
        print(ambient)

        # diffuse
        norm = glm.normalize(normal)
        light_direction = glm.normalize(light_pos - vertex_pos)
        diff = max(glm.dot(norm, light_direction), 0.0)
        diffuse = kd * diff * light_color
        print(diffuse)

        # specular
        view_dir = glm.normalize(view_pos - vertex_pos)
        reflect_dir = glm.reflect((-1) * light_direction, norm)
        spec = glm.pow(max(glm.dot(view_dir, reflect_dir), 0.0), exponent)
        specular = ks * spec * scolor
        print(specular)

        result = (ambient + diffuse + specular) * ocolor
        return result

    def drawTrianglesPhong(self, vertex_pos, indices, normals, modelT, viewT, projectionT, ocolor, scolor, k,
                           exponent,
                           lightpos, lightcolor, amb_color, doGouraud):
        i = 0
        poly_verts = []
        ocolor = glm.vec3(ocolor[0], ocolor[1], ocolor[2])
        scolor = glm.vec3(scolor[0], scolor[1], scolor[2])
        lightpos = glm.vec3(lightpos[0], lightpos[1], lightpos[2])
        lightcolor = glm.vec3(lightcolor[0], lightcolor[1], lightcolor[2])
        amb_color = glm.vec3(amb_color[0], amb_color[1], amb_color[2])
        while i < len(indices):
            # Index verts and colors according to specified stride and offsets
            index = indices[i]
            p0 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            p0_normal = glm.vec3(normals[index * 3], normals[index * 3 + 1], normals[index * 3 + 2])

            i += 1
            index = indices[i]
            p1 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            p1_normal = glm.vec3(normals[index * 3], normals[index * 3 + 1], normals[index * 3 + 2])

            i += 1
            index = indices[i]
            p2 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            p2_normal = glm.vec3(normals[index * 3], normals[index * 3 + 1], normals[index * 3 + 2])

            i += 1

            p0_normal = glm.mat3(glm.transpose(glm.inverse(modelT))) * p0_normal
            p1_normal = glm.mat3(glm.transpose(glm.inverse(modelT))) * p1_normal
            p2_normal = glm.mat3(glm.transpose(glm.inverse(modelT))) * p2_normal
            color0 = [1.0, 0.0, 1.0]
            color1 = [1.0, 0.0, 1.0]
            color2 = [1.0, 0.0, 1.0]
            if doGouraud:
                col = self.calc_lighting(modelT * p0, lightpos, p0_normal, k, ocolor, scolor, exponent,
                                         lightcolor, amb_color)
                color0 = [col[0], col[1], col[2]]
                col = self.calc_lighting(modelT * p1, lightpos, p1_normal, k, ocolor, scolor, exponent,
                                         lightcolor, amb_color)
                color1 = [col[0], col[1], col[2]]
                col = self.calc_lighting(modelT * p2, lightpos, p2_normal, k, ocolor, scolor, exponent,
                                         lightcolor, amb_color)
                color2 = [col[0], col[1], col[2]]
            screen_vert_pos0 = projectionT * viewT * modelT * p0
            screen_vert_pos1 = projectionT * viewT * modelT * p1
            screen_vert_pos2 = projectionT * viewT * modelT * p2
            vert0 = Vertex(screen_vert_pos0[0] / screen_vert_pos0[3], screen_vert_pos0[1] / screen_vert_pos0[3],
                           screen_vert_pos0[2] / screen_vert_pos0[3], None, None, color0[0], color0[1], color0[2])
            vert1 = Vertex(screen_vert_pos1[0] / screen_vert_pos1[3], screen_vert_pos1[1] / screen_vert_pos1[3],
                           screen_vert_pos1[2] / screen_vert_pos1[3], None, None, color1[0], color1[1], color1[2])
            vert2 = Vertex(screen_vert_pos2[0] / screen_vert_pos2[3], screen_vert_pos2[1] / screen_vert_pos2[3],
                           screen_vert_pos2[2] / screen_vert_pos2[3], None, None, color2[0], color2[1], color2[2])
            vert0.values["normal"] = p0_normal
            vert1.values["normal"] = p1_normal
            vert2.values["normal"] = p2_normal
            vert0.values["frags"] = modelT * p0
            vert1.values["frags"] = modelT * p1
            vert2.values["frags"] = modelT * p2
            poly = [vert0, vert1, vert2]
            print("Vertex 0 values:", vert0.values)
            print("Vertex 1 values:", vert1.values)
            print("Vertex 2 values:", vert2.values)

            # Clipping
            clipped = self.clipPoly(poly, 1.0, -1.0, 1.0, -1.0)
            poly_verts.extend(clipped)
        nverts = len(poly_verts)
        if nverts < 3:
            return
        if (nverts % 3) != 0:
            print("Bad number of vertices to define a set of triangles ", nverts, " Must be a multiple of 3")
            return
        # simply go through sets of 3 vertices and do raterization
        startidx = 0
        while startidx < nverts:
            # Add debug print statements before accessing the "normal" key later on
            print("Vertex at startidx:", poly_verts[startidx])
            print("Vertex values:", poly_verts[startidx].values)

            p0 = glm.vec3(poly_verts[startidx].x, poly_verts[startidx].y, 1)
            p0_normal = poly_verts[startidx].values["normal"]
            p0_fragpos = poly_verts[startidx].values["frags"]
            color0 = [poly_verts[startidx].r, poly_verts[startidx].g, poly_verts[startidx].b]

            sp0 = self.view_transform * p0

            poly_verts[startidx] = Vertex(int(sp0[0]), int(sp0[1]), poly_verts[startidx].z, None, None, color0[0],
                                          color0[1],
                                          color0[2])
            poly_verts[startidx].values["normal"] = p0_normal
            poly_verts[startidx].values["frags"] = p0_fragpos

            p1 = glm.vec3(poly_verts[startidx + 1].x, poly_verts[startidx + 1].y, 1)
            p1_normal = poly_verts[startidx + 1].values["normal"]
            p1_fragpos = poly_verts[startidx + 1].values["frags"]
            color1 = [poly_verts[startidx + 1].r, poly_verts[startidx + 1].g, poly_verts[startidx + 1].b]
            sp1 = self.view_transform * p1
            poly_verts[startidx + 1] = Vertex(int(sp1[0]), int(sp1[1]), poly_verts[startidx + 1].z, None, None,
                                              color1[0],
                                              color1[1], color1[2])
            poly_verts[startidx + 1].values["normal"] = p1_normal
            poly_verts[startidx + 1].values["frags"] = p1_fragpos

            p2 = glm.vec3(poly_verts[startidx + 2].x, poly_verts[startidx + 2].y, 1)
            p2_normal = poly_verts[startidx + 2].values["normal"]
            p2_fragpos = poly_verts[startidx + 2].values["frags"]
            color2 = [poly_verts[startidx + 2].r, poly_verts[startidx + 2].g, poly_verts[startidx + 2].b]
            sp2 = self.view_transform * p2
            poly_verts[startidx + 2] = Vertex(int(sp2[0]), int(sp2[1]), poly_verts[startidx + 2].z, None, None,
                                              color2[0],
                                              color2[1], color2[2])
            poly_verts[startidx + 2].values["normal"] = p2_normal
            poly_verts[startidx + 2].values["frags"] = p2_fragpos

            if doGouraud:
                self.rasterizeTriangle(poly_verts[startidx], poly_verts[startidx + 1], poly_verts[startidx + 2])
            else:
                self.rasterizeTrianglePhong(poly_verts[startidx], poly_verts[startidx + 1],
                                            poly_verts[startidx + 2],
                                            lightpos, k, ocolor, scolor, exponent, lightcolor, amb_color)
            startidx = startidx + 3

    def rasterizeTrianglePhong(self, p0, p1, p2, lightpos, k, ocolor, scolor, exponent, lightcolor, amb_color):
        min_x = min(p0.x, p1.x, p2.x)
        min_y = min(p0.y, p1.y, p2.y)
        max_x = max(p0.x, p1.x, p2.x)
        max_y = max(p0.y, p1.y, p2.y)
        # Area of triangle using cross product
        t = self.calc_cross_product(p0, p1, p2.x, p2.y)
        if t == 0:
            # Early exit in case triangle has zero area (all 3 vettices are colinear)
            return
        # Iterate over extent of triangle
        for x in range(int(min_x), int(max_x + 1)):
            for y in range(int(min_y), int(max_y + 1)):
                t2 = self.calc_cross_product(p0, p1, x, y)
                t0 = self.calc_cross_product(p1, p2, x, y)
                t1 = self.calc_cross_product(p2, p0, x, y)
                # Edge tests to verify if point lies inside the triangle (boundry inclusive)
                is_inside = (t0 >= 0) and (t1 >= 0) and (t2 >= 0)
                # Skip this pixel if it does not lie in the triangle
                if is_inside == False:
                    continue
                # Calculate barycentric coordinates
                c0 = t0 / t
                c1 = t1 / t
                c2 = t2 / t
                # interpolated fragment position
                int_x = (c2 * p2.values["frags"].x) + (c0 * p0.values["frags"].x) + (c1 * p1.values["frags"].x)
                int_y = (c2 * p2.values["frags"].y) + (c0 * p0.values["frags"].y) + (c1 * p1.values["frags"].y)
                int_z = (c2 * p2.values["frags"].z) + (c0 * p0.values["frags"].z) + (c1 * p1.values["frags"].z)
                p = glm.vec4(int_x, int_y, int_z, 1)
                # interpolated normal position
                norm_x = (c2 * p2.values["normal"].x) + (c0 * p0.values["normal"].x) + (c1 * p1.values["normal"].x)
                norm_y = (c2 * p2.values["normal"].y) + (c0 * p0.values["normal"].y) + (c1 * p1.values["normal"].y)
                norm_z = (c2 * p2.values["normal"].z) + (c0 * p0.values["normal"].z) + (c1 * p1.values["normal"].z)
                norm = glm.vec3(norm_x, norm_y, norm_z)

                col = self.calc_lighting(p, lightpos, norm, k, ocolor, scolor, exponent, lightcolor, amb_color)
                # Interpolate colors
                r = col[0]
                g = col[1]
                b = col[2]
                # Set the pixel color
                z_value = (c2 * p2.z) + (c0 * p0.z) + (c1 * p1.z)
                if z_value < self.depth[x][y]:
                    self.image[x][y] = glm.vec3(r, g, b)
                    self.depth[x][y] = z_value
                self.win.set_pixel(x, y, self.image[x][y].x, self.image[x][y].y, self.image[x][y].z)

    def drawTrianglesCheckerboard(self, vertex_pos, indices, uv, color1, color2, checksize, modelT, viewT, projectionT):
        i = 0
        poly_verts = []
        while i < len(indices):
            # Index verts and colors according to specified stride and offsets
            index = indices[i]
            p0 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            u0, v0 = uv[index * 2], uv[index * 2 + 1]
            i += 1
            index = indices[i]
            p1 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            u1, v1 = uv[index * 2], uv[index * 2 + 1]
            i += 1
            index = indices[i]
            p2 = glm.vec4(vertex_pos[index * 3], vertex_pos[index * 3 + 1], vertex_pos[index * 3 + 2], 1)
            u2, v2 = uv[index * 2], uv[index * 2 + 1]
            i += 1
            screen_vert_pos0 = projectionT * viewT * modelT * p0
            screen_vert_pos1 = projectionT * viewT * modelT * p1
            screen_vert_pos2 = projectionT * viewT * modelT * p2
            vert0 = Vertex(screen_vert_pos0[0] / screen_vert_pos0[3], screen_vert_pos0[1] / screen_vert_pos0[3],
                           screen_vert_pos0[2] / screen_vert_pos0[3], u0, v0, None, None, None)
            vert1 = Vertex(screen_vert_pos1[0] / screen_vert_pos1[3], screen_vert_pos1[1] / screen_vert_pos1[3],
                           screen_vert_pos1[2] / screen_vert_pos1[3], u1, v1, None, None, None)
            vert2 = Vertex(screen_vert_pos2[0] / screen_vert_pos2[3], screen_vert_pos2[1] / screen_vert_pos2[3],
                           screen_vert_pos2[2] / screen_vert_pos2[3], u2, v2, None, None, None)
            poly = [vert0, vert1, vert2]
            # Clipping
            clipped = self.clipPoly(poly, 1.0, -1.0, 1.0, -1.0)
            poly_verts.extend(clipped)
        nverts = len(poly_verts)
        if nverts < 3:
            return
        if (nverts % 3) != 0:
            print("Bad number of verticies to define a set of triangles ", nverts, " Must be a multiple of 3")
            return
        # simply go through sets of 3 vertices and do raterization
        startidx = 0
        while startidx < nverts:
            p0 = glm.vec3(poly_verts[startidx].x, poly_verts[startidx].y, 1)
            # color0 = [poly_verts[startidx].r, poly_verts[startidx].g, poly_verts[startidx].b]
            uv0 = [poly_verts[startidx].u, poly_verts[startidx].v]
            sp0 = self.view_transform * p0
            poly_verts[startidx] = Vertex(int(sp0[0]), int(sp0[1]), poly_verts[startidx].z, uv0[0], uv0[1], None, None,
                                          None)

            p1 = glm.vec3(poly_verts[startidx + 1].x, poly_verts[startidx + 1].y, 1)
            # color1 = [poly_verts[startidx + 1].r, poly_verts[startidx + 1].g, poly_verts[startidx + 1].b]
            uv1 = [poly_verts[startidx + 1].u, poly_verts[startidx + 1].v]
            sp1 = self.view_transform * p1
            poly_verts[startidx + 1] = Vertex(int(sp1[0]), int(sp1[1]), poly_verts[startidx + 1].z, uv1[0],
                                              uv1[1], None, None, None)

            p2 = glm.vec3(poly_verts[startidx + 2].x, poly_verts[startidx + 2].y, 1)
            # color2 = [poly_verts[startidx + 2].r, poly_verts[startidx + 2].g, poly_verts[startidx + 2].b]
            uv2 = [poly_verts[startidx + 2].u, poly_verts[startidx + 2].v]
            sp2 = self.view_transform * p2
            poly_verts[startidx + 2] = Vertex(int(sp2[0]), int(sp2[1]), poly_verts[startidx + 2].z, uv2[0],
                                              uv2[1], None, None, None)
            self.rasterizeTriangleCheckered(poly_verts[startidx], poly_verts[startidx + 1], poly_verts[startidx + 2],
                                            color1, color2, checksize)
            startidx = startidx + 3

    def rasterizeTriangleCheckered(self, p0, p1, p2, color1, color2, checksize):
        min_x = min(p0.x, p1.x, p2.x)
        min_y = min(p0.y, p1.y, p2.y)
        max_x = max(p0.x, p1.x, p2.x)
        max_y = max(p0.y, p1.y, p2.y)
        # Area of triangle using cross product
        t = self.calc_cross_product(p0, p1, p2.x, p2.y)
        if t == 0:
            # Early exit in case triangle has zero area (all 3 vettices are colinear)
            return
        # Iterate over extent of triangle
        for x in range(int(min_x), int(max_x + 1)):
            for y in range(int(min_y), int(max_y + 1)):
                t2 = self.calc_cross_product(p0, p1, x, y)
                t0 = self.calc_cross_product(p1, p2, x, y)
                t1 = self.calc_cross_product(p2, p0, x, y)
                # Edge tests to verify if point lies inside the triangle (boundry inclusive)
                is_inside = (t0 >= 0) and (t1 >= 0) and (t2 >= 0)
                # Skip this pixel if it does not lie in the triangle
                if is_inside == False:
                    continue
                # Calculate barycentric coordinates
                c0 = t0 / t
                c1 = t1 / t
                c2 = t2 / t
                # Interpolate colors
                u = (c0 * p0.u) + (c1 * p1.u) + (c2 * p2.u)
                v = (c0 * p0.v) + (c1 * p1.v) + (c2 * p2.v)
                print("Pixel at (", x, ",", y, ") - UV coordinates: (", u, ",", v, ")")
                urow = int(u / checksize)
                vrow = int(v / checksize)

                if (urow % 2 == 0 and vrow % 2 == 0) or (urow % 2 != 0 and vrow % 2 != 0):
                    c_r = color1[0]
                    c_g = color1[1]
                    c_b = color1[2]
                else:
                    c_r = color2[0]
                    c_g = color2[1]
                    c_b = color2[2]

                # Set the pixel color
                z_value = (c2 * p2.z) + (c0 * p0.z) + (c1 * p1.z)
                if z_value < self.depth[x][y]:
                    self.image[x][y] = glm.vec3(c_r, c_g, c_b)
                    self.depth[x][y] = z_value
                self.win.set_pixel(x, y, self.image[x][y].x, self.image[x][y].y, self.image[x][y].z)
