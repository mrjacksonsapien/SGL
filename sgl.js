/**
 * This is a basic, barely-optimized 3D engine where code clarity has been prioritized over performance.
 * The goal of this library is to fully illustrate in a more human-readable manner all the steps of 3D rendering
 * pipeline by using only sequential instructions (CPU-based instructions) which means sacrificing the parallelism
 * of the GPU. Warning: This library is made for learning purpose of the 3D graphic pipeline and not for heavy intensive
 * scenes rendering. The CanvasRenderingContext2D API is used for the final render.
 *
 * Important note: You can use different data structures but most of the time, functions assume that you pass as an
 * argument a Float32Array and will directly modify the elements inside it.
 *
 * @author mrjacksonsapien
 */

/**
 * Math functions.
 */
export class SGLMath {
    /**
     * Converts degrees to radians.
     * @param degrees Unit in degrees.
     * @returns {number} Unit in radians.
     */
    static degToRad(degrees) {
        return degrees * (Math.PI / 180);
    }

    /**
     * Cotangent of opening of frustum.
     * @param degrees The field of view of the camera frustum
     * @returns {number} Cotangent
     */
    static cot(degrees) {
        return 1 / Math.tan(SGLMath.degToRad(degrees));
    }
}

/**
 * Takes care of rendering. Uses the Matrix class for vertices transformation.
 */
export class Renderer {
    /**
     * The size of one triangle inside the triangle data array.
     * @type {number}
     */
    static SIZEOF_TRIANGLE_DATA = 6;
    /**
     * The size of one vertex inside the vertex data array.
     * @type {number}
     */
    static SIZEOF_VERTEX_DATA = 4;

    /**
     * Creates an instance to render the scene on the canvas using CanvasRenderingContext2D API.
     * @param canvas HTML canvas.
     * @param scene SGL Scene object.
     */
    constructor(canvas, scene) {
        this.canvas = canvas;
        this.scene = scene;
    }

    /**
     * Converts the Scene object attached to the Renderer into 2 Float32Arrays stored in one JavaScript array.
     * @returns {Float32Array[]} Contains triangles data and vertices data. The triangles array contains the indices
     * to their corresponding vertices in the vertices data array.
     */
    convertSceneToFlatArrays() {
        let vertices = [];
        const triangles = [];
        const meshes = this.scene.meshes;

        // Get all triangles and vertices from the instances
        for (let i = 0; i < meshes.length; i++) {
            const mesh = meshes[i];

            for (let j = 0; j < mesh.triangles.length; j++) {
                triangles.push(mesh.triangles[j]);
            }

            for (let j = 0; j < mesh.vertices.length; j++) {
                vertices.push(mesh.vertices[j]);
            }
        }

        const trianglesData = new Float32Array(triangles.length * Renderer.SIZEOF_TRIANGLE_DATA);

        // Map vertices to triangles
        for (let i = 0; i < triangles.length; i++) {
            const triangle = triangles[i];
            const triangleIndex = i * Renderer.SIZEOF_TRIANGLE_DATA;

            trianglesData[triangleIndex] = vertices.indexOf(triangle.vertex1) * Renderer.SIZEOF_VERTEX_DATA;
            trianglesData[triangleIndex + 1] = vertices.indexOf(triangle.vertex2) * Renderer.SIZEOF_VERTEX_DATA;
            trianglesData[triangleIndex + 2] = vertices.indexOf(triangle.vertex3) * Renderer.SIZEOF_VERTEX_DATA;
            trianglesData[triangleIndex + 3] = triangle.color.r;
            trianglesData[triangleIndex + 4] = triangle.color.g;
            trianglesData[triangleIndex + 5] = triangle.color.b;
        }

        let verticesData = new Float32Array(vertices.length * Renderer.SIZEOF_VERTEX_DATA);

        // Add vertex values into flat array
        for (let i = 0; i < vertices.length; i++) {
            const vertex = vertices[i];
            const verticesDataIndex = i * Renderer.SIZEOF_VERTEX_DATA;

            verticesData[verticesDataIndex] = vertex.position.x;
            verticesData[verticesDataIndex + 1] = vertex.position.y;
            verticesData[verticesDataIndex + 2] = vertex.position.z;
            verticesData[verticesDataIndex + 3] = 1;
        }

        return [trianglesData, verticesData];
    }

    /**
     * Apply a matrix to a flat array of many vertices (Converts each vertex individually).
     * @param matrix The matrix obtained from a method in the Matrix class.
     * @param verticesData The vertices data float array.
     */
    applyMatrixToVertices(matrix, verticesData) {
        for (let i = 0; i < verticesData.length / Renderer.SIZEOF_VERTEX_DATA; i++) {
            Matrix.multiplyMatrixWithVertex(matrix,  i * Renderer.SIZEOF_VERTEX_DATA, verticesData);
        }
    }

    /**
     * Culls triangles that are facing away from camera (discards them). Front or back of the triangle is defined
     * by the order the vertices are defined in the triangle.
     * @param trianglesData
     * @param verticesData
     * @returns {Float32Array[]}
     */
    cull(trianglesData, verticesData) {
        const verticesIndexMap = {};
        const keptTriangles = [];
        const keptVertices = [];

        function addVertex(verticesDataIndex) {
            if (verticesIndexMap[verticesDataIndex] === undefined) {
                const nextIndex = keptVertices.length;

                keptVertices.push(
                    verticesData[verticesDataIndex],
                    verticesData[verticesDataIndex + 1],
                    verticesData[verticesDataIndex + 2],
                    verticesData[verticesDataIndex + 3]
                )

                verticesIndexMap[verticesDataIndex] = nextIndex;
            }
            return verticesIndexMap[verticesDataIndex];
        }

        function getXYZ(vertexIndex) {
            return [
                verticesData[vertexIndex],
                verticesData[vertexIndex + 1],
                verticesData[vertexIndex + 2]
            ];
        }

        for (let i = 0; i < trianglesData.length / Renderer.SIZEOF_TRIANGLE_DATA; i++) {
            const triangleIndex = i * Renderer.SIZEOF_TRIANGLE_DATA;

            const vertex1Index = trianglesData[triangleIndex];
            const vertex2Index = trianglesData[triangleIndex + 1];
            const vertex3Index = trianglesData[triangleIndex + 2];

            const vertex1XYZ = getXYZ(vertex1Index);
            const vertex2XYZ = getXYZ(vertex2Index);
            const vertex3XYZ = getXYZ(vertex3Index);

            // Counter clock-wise convention (CCW)
            const edgeAX = vertex2XYZ[0] - vertex1XYZ[0];
            const edgeAY = vertex2XYZ[1] - vertex1XYZ[1];
            const edgeAZ = vertex2XYZ[2] - vertex1XYZ[2];

            const edgeBX = vertex3XYZ[0] - vertex1XYZ[0];
            const edgeBY = vertex3XYZ[1] - vertex1XYZ[1];
            const edgeBZ = vertex3XYZ[2] - vertex1XYZ[2];

            const normalX = (edgeAY * edgeBZ) - (edgeAZ * edgeBY);
            const normalY = (edgeAZ * edgeBX) - (edgeAX * edgeBZ);
            const normalZ = (edgeAX * edgeBY) - (edgeAY * edgeBX);

            const viewDirectionX = vertex1XYZ[0];
            const viewDirectionY = vertex1XYZ[1];
            const viewDirectionZ = vertex1XYZ[2];

            const dot = (normalX * viewDirectionX) + (normalY * viewDirectionY) + (normalZ * viewDirectionZ);

            if (dot < 0) {
                keptTriangles.push(
                    addVertex(vertex1Index),
                    addVertex(vertex2Index),
                    addVertex(vertex3Index),
                    trianglesData[triangleIndex + 3],
                    trianglesData[triangleIndex + 4],
                    trianglesData[triangleIndex + 5]
                )
            }
        }

        return [new Float32Array(keptTriangles), new Float32Array(keptVertices)];
    }

    /**
     * Clips all the triangles and returns a new JavaScript array containing new data for the
     * triangles and vertices (same structure as return of convertSceneToFlatArrays()).
     * @param trianglesData The triangles data.
     * @param verticesData The vertices data.
     * @returns {Float32Array[]} The new JavaScript array containing the triangles and vertices data.
     */
    clip(trianglesData, verticesData) {
        let activeTriangles = [...trianglesData];
        let activeVertices = [...verticesData];

        let nextTriangles = [];
        let nextVertices = [];

        function planesRelation(vertexIndex, i) {
            const x = activeVertices[vertexIndex];
            const y = activeVertices[vertexIndex + 1];
            const z = activeVertices[vertexIndex + 2];
            const w = activeVertices[vertexIndex + 3];

            const planes = [
                x >= w, // Left
                x <= -w, // Right
                y >= w, // Top
                y <= -w, // Bottom
                z >= w, // Near
                z <= -w // Far
            ]

            return planes[i];
        }

        function getIntersection(p1Index, p2Index, planeIndex) {
            // The w coefficient is inverted because the matrices respect the OpenGL standard.
            const planesConstant = [
                [1, 0, 0, -1], // Left
                [-1, 0, 0, -1], // Right
                [0, 1, 0, -1], // Bottom
                [0, -1, 0, -1], // Top
                [0, 0, 1, -1], // Near
                [0, 0, -1, -1] // Far
            ]

            const plane = planesConstant[planeIndex];

            const x1 = activeVertices[p1Index];
            const y1 = activeVertices[p1Index + 1];
            const z1 = activeVertices[p1Index + 2];
            const w1 = activeVertices[p1Index + 3];

            const x2 = activeVertices[p2Index];
            const y2 = activeVertices[p2Index + 1];
            const z2 = activeVertices[p2Index + 2];
            const w2 = activeVertices[p2Index + 3];

            const numerator = -(plane[0] * x1 + plane[1] * y1 + plane[2] * z1 + plane[3] * w1);
            const denominator = (plane[0] * x2 + plane[1] * y2 + plane[2] * z2 + plane[3] * w2) + numerator;

            if (Math.abs(denominator) < 1e-6) {
                return null;
            }

            const t = numerator / denominator;

            const x = x1 + t * (x2 - x1);
            const y = y1 + t * (y2 - y1);
            const z = z1 + t * (z2 - z1);
            const w = w1 + t * (w2 - w1);

            return [x, y, z, w];
        }

        // Iterate over each plane
        for (let i = 0; i < 6; i++) {
            let verticesIndexMap = {};

            function addVertex(activeIndex) {
                if (verticesIndexMap[activeIndex] === undefined) {
                    const nextIndex = nextVertices.length;

                    nextVertices.push(
                        activeVertices[activeIndex],
                        activeVertices[activeIndex + 1],
                        activeVertices[activeIndex + 2],
                        activeVertices[activeIndex + 3]
                    )

                    verticesIndexMap[activeIndex] = nextIndex;
                }
                return verticesIndexMap[activeIndex];
            }

            function createVertex(vertexData) {
                if (vertexData.length === Renderer.SIZEOF_VERTEX_DATA) {
                    const nextIndex = nextVertices.length;
                    nextVertices.push(...vertexData);
                    return nextIndex;
                }
            }

            function createNewTriangle(vertex1Index, vertex2Index, vertex3Index, originalTriangleIndex) {
                const newTriangleIndex = nextTriangles.length;
                nextTriangles.push(
                    vertex1Index,
                    vertex2Index,
                    vertex3Index,
                    activeTriangles[originalTriangleIndex + 3],
                    activeTriangles[originalTriangleIndex + 4],
                    activeTriangles[originalTriangleIndex + 5]
                )
                return newTriangleIndex;
            }

            for (let j = 0; j < activeTriangles.length / Renderer.SIZEOF_TRIANGLE_DATA; j++) {
                const triangleIndex = j * Renderer.SIZEOF_TRIANGLE_DATA;
                const vertex1Index = activeTriangles[triangleIndex];
                const vertex2Index = activeTriangles[triangleIndex + 1];
                const vertex3Index = activeTriangles[triangleIndex + 2];

                const inside = [];
                const outside = [];

                if (planesRelation(vertex1Index, i)) inside.push(vertex1Index)
                else outside.push(vertex1Index)

                if (planesRelation(vertex2Index, i)) inside.push(vertex2Index)
                else outside.push(vertex2Index)

                if (planesRelation(vertex3Index, i)) inside.push(vertex3Index)
                else outside.push(vertex3Index)

                if (inside.length === 3) {
                    createNewTriangle(
                        addVertex(inside[0]),
                        addVertex(inside[1]),
                        addVertex(inside[2]),
                        triangleIndex
                    );
                } else if (inside.length === 2) {
                    const intersection1 = createVertex(getIntersection(inside[0], outside[0], i));
                    const intersection2 = createVertex(getIntersection(inside[1], outside[0], i));

                    createNewTriangle(
                        addVertex(inside[0]),
                        addVertex(inside[1]),
                        intersection1,
                        triangleIndex
                    );

                    createNewTriangle(
                        addVertex(inside[1]),
                        intersection2,
                        intersection1,
                        triangleIndex
                    );
                } else if (inside.length === 1) {
                    const intersection1 = createVertex(getIntersection(inside[0], outside[0], i));
                    const intersection2 = createVertex(getIntersection(inside[0], outside[1], i));

                    createNewTriangle(
                        addVertex(inside[0]),
                        intersection2,
                        intersection1,
                        triangleIndex
                    );
                }
            }

            activeTriangles = nextTriangles;
            activeVertices = nextVertices;
            nextTriangles = [];
            nextVertices = [];
        }

        return [new Float32Array(activeTriangles), new Float32Array(activeVertices)];
    }

    /**
     * Applies perspective division to clip-space vertices.
     * @param verticesData Clip-space vertices.
     */
    applyPerspectiveDivisionToClipVertices(verticesData) {
        for (let i = 0; i < verticesData.length / Renderer.SIZEOF_VERTEX_DATA; i++) {
            const vertexIndex = i * Renderer.SIZEOF_VERTEX_DATA;

            verticesData[vertexIndex] /= verticesData[vertexIndex + 3];
            verticesData[vertexIndex + 1] /= verticesData[vertexIndex + 3];
            verticesData[vertexIndex + 2] /= verticesData[vertexIndex + 3];
        }
    }

    /**
     * Maps the NDC vertices to the canvas 2D coordinates.
     * @param verticesData NDC vertices.
     */
    mapNdcVerticesToScreenCoordinates(verticesData) {
        for (let i = 0; i < verticesData.length / Renderer.SIZEOF_VERTEX_DATA; i++) {
            const vertexIndex = i * Renderer.SIZEOF_VERTEX_DATA;

            verticesData[vertexIndex] = (verticesData[vertexIndex] + 1) / 2 * this.canvas.clientWidth;
            verticesData[vertexIndex + 1] = (1 - verticesData[vertexIndex + 1]) / 2 * this.canvas.clientHeight;
        }
    }

    renderTriangles(trianglesData, verticesData, ctx) {
        function pointIsInTriangle(px, py, ax, ay, bx, by, cx, cy, isCCW) {
            let edge1 = (py - ay) * (bx - ax) - (px - ax) * (by - ay); // Edge AB
            let edge2 = (py - by) * (cx - bx) - (px - bx) * (cy - by); // Edge BC
            let edge3 = (py - cy) * (ax - cx) - (px - cx) * (ay - cy); // Edge CA

            if (!isCCW) {
                edge1 *= -1;
                edge2 *= -1;
                edge3 *= -1;
            }

            return (edge1 >= 0) && (edge2 >= 0) && (edge3 >= 0);
        }

        function isTriangleCCW(ax, ay, bx, by, cx, cy) {
            // Calculate the signed area
            const signedArea = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
            return signedArea > 0; // CCW if positive
        }

        ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        const imageData = ctx.getImageData(0, 0, this.canvas.width, this.canvas.height);
        const depthBuffer = new Float32Array(this.canvas.width * this.canvas.height);
        depthBuffer.fill(Infinity);

        for (let i = 0; i < trianglesData.length / Renderer.SIZEOF_TRIANGLE_DATA; i++) {
            const triangleIndex = i * Renderer.SIZEOF_TRIANGLE_DATA;

            const vertex1Index = trianglesData[triangleIndex];
            const vertex2Index = trianglesData[triangleIndex + 1];
            const vertex3Index = trianglesData[triangleIndex + 2];

            const vertex1X = verticesData[vertex1Index];
            const vertex1Y = verticesData[vertex1Index + 1];

            const vertex2X = verticesData[vertex2Index];
            const vertex2Y = verticesData[vertex2Index + 1];

            const vertex3X = verticesData[vertex3Index];
            const vertex3Y = verticesData[vertex3Index + 1];

            const isCCW = isTriangleCCW(vertex1X, vertex1Y, vertex2X, vertex2Y, vertex3X, vertex3Y);

            const minX = Math.floor(Math.min(vertex1X, vertex2X, vertex3X));
            const maxX = Math.ceil(Math.max(vertex1X, vertex2X, vertex3X));
            const minY = Math.floor(Math.min(vertex1Y, vertex2Y, vertex3Y));
            const maxY = Math.ceil(Math.max(vertex1Y, vertex2Y, vertex3Y));

            for (let x = minX; x < maxX; x++) {
                for (let y = minY; y < maxY; y++) {
                    if (pointIsInTriangle(x, y, vertex1X, vertex1Y, vertex2X, vertex2Y, vertex3X, vertex3Y, isCCW)) {
                        const vertex1Z = verticesData[vertex1Index + 2];
                        const vertex2Z = verticesData[vertex2Index + 2];
                        const vertex3Z = verticesData[vertex3Index + 2];

                        // Compute barycentric coordinates (α1, α2, α3)
                        const denominator = (vertex2Y - vertex3Y) * (vertex1X - vertex3X) + (vertex3X - vertex2X) * (vertex1Y - vertex3Y);

                        const a1 = ((vertex2Y - vertex3Y) * (x - vertex3X) + (vertex3X - vertex2X) * (y - vertex3Y)) / denominator;
                        const a2 = ((vertex3Y - vertex1Y) * (x - vertex3X) + (vertex1X - vertex3X) * (y - vertex3Y)) / denominator;
                        const a3 = 1 - a1 - a2;

                        // Interpolate z value using barycentric coordinates
                        const z = a1 * vertex1Z + a2 * vertex2Z + a3 * vertex3Z;

                        const pixelIndex = (y * this.canvas.width + x);
                        if (z < depthBuffer[pixelIndex]) {
                            depthBuffer[pixelIndex] = z;
                            const data = imageData.data;
                            const index = (y * this.canvas.width + x) * 4;

                            data[index] = trianglesData[triangleIndex + 3];
                            data[index + 1] = trianglesData[triangleIndex + 4];
                            data[index + 2] = trianglesData[triangleIndex + 5];
                            data[index + 3] = 255;
                        }
                    }
                }
            }
        }

        ctx.putImageData(imageData, 0, 0);
    }

    /**
     * Goes through all the steps of the 3D rendering pipeline. The scene must have a current camera when this method is called.
     * Since JavaScript is really easy to break, make sure that all your attributes are initialized and defined. If you use the classes
     * correctly, it shouldn't be a problem but be careful.
     */
    render() {
        const camera = this.scene.currentCamera;

        let sceneData = this.convertSceneToFlatArrays();

        // Now in world space

        this.applyMatrixToVertices(Matrix.createViewMatrix(camera), sceneData[1]);

        // Now in view space
        sceneData = this.cull(sceneData[0], sceneData[1]);

        this.applyMatrixToVertices(Matrix.createProjectionMatrix(camera), sceneData[1]);

        // Now in clip space

        sceneData = this.clip(sceneData[0], sceneData[1]);

        this.applyPerspectiveDivisionToClipVertices(sceneData[1]);

        // Now in NDC space

        this.mapNdcVerticesToScreenCoordinates(sceneData[1]);

        // Now in 2D screen space
        const ctx = this.canvas.getContext('2d', {willReadFrequently: true});

        this.renderTriangles(sceneData[0], sceneData[1], ctx);
    }
}

/**
 * Holds all the hardcoded matrices for 3D transformation. Some conventional aspects are respected while others like
 * matrix multiplications might change from one engine to another.
 */
export class Matrix {
    /**
     * Multiplies a 4x4 matrix (flat array) with a vertex.
     * @param m The matrix (flat array).
     * @param vertexIndex The index of the first component of the vertex (x) inside the verticesData.
     * @param verticesData The vertices data (Float32Array).
     */
    static multiplyMatrixWithVertex(m, vertexIndex, verticesData) {
        const x = verticesData[vertexIndex];
        const y = verticesData[vertexIndex + 1];
        const z = verticesData[vertexIndex + 2];
        const w = verticesData[vertexIndex + 3];

        verticesData[vertexIndex] = x * m[0] + y * m[4] + z * m[8] + w * m[12];
        verticesData[vertexIndex + 1] = x * m[1] + y * m[5] + z * m[9] + w * m[13];
        verticesData[vertexIndex + 2] = x * m[2] + y * m[6] + z * m[10] + w * m[14];
        verticesData[vertexIndex + 3] = x * m[3] + y * m[7] + z * m[11] + w * m[15];
    }

    /**
     * Multiplies two 4x4 matrices together.
     * @param a First matrix (flat array).
     * @param b Second matrix (flat array).
     * @returns {Float32Array} Result 4x4 matrix.
     */
    static multiply4x4Matrices(a, b) {
        const result = new Float32Array(16);

        for (let i = 0; i < 4; i++) {
            const offsetA = i * 4;
            for (let j = 0; j < 4; j++) {
                const offsetB = j;
                result[offsetA + j] = a[offsetA] * b[offsetB] +
                    a[offsetA + 1] * b[4 + offsetB] +
                    a[offsetA + 2] * b[8 + offsetB] +
                    a[offsetA + 3] * b[12 + offsetB];
            }
        }

        return result;
    }

    /**
     * Creates a translation matrix with x, y, z position components.
     * @param x
     * @param y
     * @param z
     * @returns {Float32Array}
     */
    static createTranslationMatrix(x, y, z) {
        const translationMatrix = new Float32Array(16);
        translationMatrix.set([
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            x, y, z, 1
        ]);
        return translationMatrix;
    }

    /**
     * Creates a pitch matrix (x rotation).
     * @param pitchDegrees Degrees for x.
     * @returns {Float32Array}
     */
    static createPitchMatrix(pitchDegrees) {
        const pitch = SGLMath.degToRad(pitchDegrees);
        const pitchMatrix = new Float32Array(16);
        pitchMatrix.set([
            1, 0, 0, 0,
            0, Math.cos(pitch), -Math.sin(pitch), 0,
            0, Math.sin(pitch), Math.cos(pitch), 0,
            0, 0, 0, 1
        ]);
        return pitchMatrix;
    }

    /**
     * Create a yaw matrix (y).
     * @param yawDegrees Degrees for y.
     * @returns {Float32Array}
     */
    static createYawMatrix(yawDegrees) {
        const yaw = SGLMath.degToRad(yawDegrees);
        const yawMatrix = new Float32Array(16);
        yawMatrix.set([
            Math.cos(yaw), 0, Math.sin(yaw), 0,
            0, 1, 0, 0,
            -Math.sin(yaw), 0, Math.cos(yaw), 0,
            0, 0, 0, 1
        ]);
        return yawMatrix;
    }

    /**
     * Creates a roll matrix (z).
     * @param rollDegrees Degrees for z.
     * @returns {Float32Array}
     */
    static createRollMatrix(rollDegrees) {
        const roll = SGLMath.degToRad(rollDegrees);
        const rollMatrix = new Float32Array(16);
        rollMatrix.set([
            Math.cos(roll), -Math.sin(roll), 0, 0,
            Math.sin(roll), Math.cos(roll), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        ]);
        return rollMatrix;
    }

    /**
     * Creates an euleur matrix respecting the video game camera rotation standard.
     * @param x In degrees.
     * @param y In degrees.
     * @param z In degrees.
     * @returns {Float32Array}
     */
    static createEulerMatrix(x, y, z) {
        return Matrix.multiply4x4Matrices(
            Matrix.multiply4x4Matrices(
                Matrix.createYawMatrix(y),
                Matrix.createPitchMatrix(x)
            ),
            Matrix.createRollMatrix(z)
        );
    }

    /**
     * Creates a view matrix for the camera.
     * @param camera Camera object.
     * @returns {Float32Array}
     */
    static createViewMatrix(camera) {
        return Matrix.multiply4x4Matrices(
            Matrix.createTranslationMatrix(camera.position.x, camera.position.y, camera.position.z),
            Matrix.createEulerMatrix(camera.orientation.x, camera.orientation.y, camera.orientation.z)
        )
    }

    /**
     * Creates a projection matrix (clip-space transformation).
     * @param camera Camera object.
     * @returns {Float32Array}
     */
    static createProjectionMatrix(camera) {
        const aspectRatio = camera.canvas.clientWidth / camera.canvas.clientHeight;
        const projectionMatrix = new Float32Array(16);
        projectionMatrix.set([
            SGLMath.cot(camera.fov / 2) / aspectRatio, 0, 0, 0,
            0, SGLMath.cot(camera.fov / 2), 0, 0,
            0, 0, -(camera.far / (camera.far - camera.near)), -1,
            0, 0, (camera.far * camera.near) / (camera.far - camera.near), 0
        ]);
        return projectionMatrix;
    }
}

/**
 * Contains all the meshes that will be rendered and the current camera that will be used to render the scene.
 */
export class Scene {
    /**
     * The camera that will be set as the current camera for rendering.
     * @param camera Camera object.
     */
    constructor(camera) {
        this.currentCamera = camera;
        this.meshes = [];
    }

    /**
     * Add a mesh to the scene.
     * @param mesh
     */
    add(mesh) {
        this.meshes.push(mesh);
    }

    /**
     * Removes a mesh from the scene.
     * @param mesh
     */
    remove(mesh) {
        let index = this.meshes.indexOf(mesh);
        if (index !== -1) {
            this.meshes.splice(index, 1);
        }
    }

    /**
     * Logic to update the scene.
     */
    act() {
        this.meshes.forEach(mesh => {
            mesh.act();
        });
    }
}

/**
 * Camera object used for rendering a scene.
 */
export class Camera {
    /**
     * Creates the camera.
     * @param canvas HTML canvas used for rendering the scene.
     * @param near The near clipping plane (cannot be equal or under 0).
     * @param far Has to be larger than the near value.
     * @param fov Field of view.
     * @param position Position of the camera.
     * @param orientation Orientation of the camera.
     */
    constructor(canvas, near, far, fov, position, orientation) {
        this.canvas = canvas;
        this.near = near;
        this.far = far;
        this.fov = fov;
        this.position = position;
        this.orientation = orientation;
    }
}

/**
 * Mesh object containing all the data for rendering.
 */
export class Mesh {
    constructor(vertices, triangles) {
        this.vertices = vertices;
        this.triangles = triangles;
        this.act = () => {};
    }
}

/**
 * Contains the references to 3 Vertex and the data of the triangle.
 */
export class Triangle {
    constructor(vertex1, vertex2, vertex3, color) {
        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
        this.vertex3 = vertex3;
        this.color = color;
    }
}

/**
 * Contains the data of a vertex.
 */
export class Vertex {
    constructor(position) {
        this.position = position;
    }
}

/**
 * Contains an x, y and z component.
 */
export class Vector3 {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

export class Vector2 {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }
}

export class Color {
    static RED = new Color(255, 0, 0);
    static GREEN = new Color(0, 255, 0);
    static BLUE = new Color(0, 0, 255);

    constructor(r, g, b) {
        this.r = r;
        this.g = g;
        this.b = b;
    }
}

// Templates //

export class Cube extends Mesh {
    constructor(position, size) {
        super();
        this.position = position;
        this.size = size;

        this.vertices = [
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y - this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y - this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y + this.size.y / 2, this.position.z + this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y + this.size.y / 2, this.position.z + this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y - this.size.y / 2, this.position.z + this.size.z / 2)),
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y - this.size.y / 2, this.position.z + this.size.z / 2))
        ];

        this.triangles = [
            new Triangle(this.vertices[0], this.vertices[1], this.vertices[2], Color.BLUE),
            new Triangle(this.vertices[0], this.vertices[2], this.vertices[3], Color.BLUE),
            new Triangle(this.vertices[6], this.vertices[5], this.vertices[4], Color.BLUE),
            new Triangle(this.vertices[7], this.vertices[6], this.vertices[4], Color.BLUE),
            new Triangle(this.vertices[5], this.vertices[1], this.vertices[0], Color.RED),
            new Triangle(this.vertices[4], this.vertices[5], this.vertices[0], Color.RED),
            new Triangle(this.vertices[3], this.vertices[2], this.vertices[6], Color.RED),
            new Triangle(this.vertices[3], this.vertices[6], this.vertices[7], Color.RED),
            new Triangle(this.vertices[3], this.vertices[4], this.vertices[0], Color.GREEN),
            new Triangle(this.vertices[4], this.vertices[3], this.vertices[7], Color.GREEN),
            new Triangle(this.vertices[1], this.vertices[5], this.vertices[6], Color.GREEN),
            new Triangle(this.vertices[6], this.vertices[2], this.vertices[1], Color.GREEN)
        ];
    }
}

export class TriangleMesh extends Mesh {
    constructor(position, size) {
        super();
        this.position = position;
        this.size = size;

        this.vertices = [
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y - this.size.y / 2 + 1, this.position.z - this.size.z / 2 - 1))
        ];

        this.triangles = [
            new Triangle(this.vertices[0], this.vertices[1], this.vertices[2], Color.RED)
        ];
    }
}