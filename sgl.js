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

        // Get all the vertices from all the triangles in the scene
        for (let i = 0; i < meshes.length; i++) {
            const mesh = meshes[i];

            for (let j = 0; j < mesh.triangles.length; j++) {
                const triangle = mesh.triangles[j];
                triangles.push(triangle);

                vertices.push(triangle.vertex1);
                vertices.push(triangle.vertex2);
                vertices.push(triangle.vertex3);
            }
        }

        // Remove duplicates (many triangles sharing the same vertex)
        vertices = [...new Set(vertices)];

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
                nextTriangles.push(
                    vertex1Index,
                    vertex2Index,
                    vertex3Index,
                    activeTriangles[originalTriangleIndex + 3],
                    activeTriangles[originalTriangleIndex + 4],
                    activeTriangles[originalTriangleIndex + 5]
                )
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
                        intersection1,
                        intersection2,
                        triangleIndex
                    );
                } else if (inside.length === 1) {
                    const intersection1 = createVertex(getIntersection(inside[0], outside[0], i));
                    const intersection2 = createVertex(getIntersection(inside[0], outside[1], i));

                    createNewTriangle(
                        addVertex(inside[0]),
                        intersection1,
                        intersection2,
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

    /**
     * Display the triangle's wireframe.
     * @param trianglesData Triangles data after clipping.
     * @param verticesData Vertices data after 2D mapping.
     */
    displayTrianglesWireframe(trianglesData, verticesData) {
        const ctx = this.canvas.getContext('2d');
        for (let i = 0; i < trianglesData.length / Renderer.SIZEOF_TRIANGLE_DATA; i++) {
            const triangleIndex = i * Renderer.SIZEOF_TRIANGLE_DATA;

            const vertex1Index = trianglesData[triangleIndex];
            const vertex2Index = trianglesData[triangleIndex + 1];
            const vertex3Index = trianglesData[triangleIndex + 2];

            ctx.beginPath();
            ctx.moveTo(verticesData[vertex1Index], verticesData[vertex1Index + 1]);
            ctx.lineTo(verticesData[vertex2Index], verticesData[vertex2Index + 1]);
            ctx.lineTo(verticesData[vertex3Index], verticesData[vertex3Index + 1]);
            ctx.closePath();
            ctx.fillStyle = `rgb(${trianglesData[triangleIndex + 3]}, ${trianglesData[triangleIndex + 4]}, ${trianglesData[triangleIndex + 5]})`;
            ctx.fill();
        }
    }

    /**
     * Goes through all the steps of the 3D rendering pipeline. The scene must have a current camera when this method is called.
     * Since JavaScript is really easy to break, make sure that all your attributes are initialized and defined. If you use the classes
     * correctly, it shouldn't be a problem but be careful.
     */
    render() {
        const camera = this.scene.currentCamera;

        let sceneData = this.convertSceneToFlatArrays();
        let trianglesData = sceneData[0];
        let verticesData = sceneData[1];

        /* Now in world space (coordinates relative to the world origin). The scene has been converted into
        2 flat (one dimension) arrays: One for all the vertices in the scene, and one for the triangles. The triangles
        hold the indices (in verticesData) of the vertices that shapes them. You can see this as some kind of "reference".
        Making flat arrays will help make calculations faster than accessing objects data. */

        this.applyMatrixToVertices(Matrix.createViewMatrix(camera), verticesData);

        // Now in view space: Perform calculations for lighting, normal vectors, etc.

        this.applyMatrixToVertices(Matrix.createProjectionMatrix(camera), verticesData);

        /* Now in clip space: Perform clipping and triangles subdivision. The variables holding the data arrays are updated
        because some triangles and vertices might have been discarded with some new ones added, which requires changing the size
        of the arrays. */

        sceneData = this.clip(trianglesData, verticesData);
        trianglesData = sceneData[0];
        verticesData = sceneData[1];

        this.applyPerspectiveDivisionToClipVertices(verticesData);

        /* Now in NDC space. The perspective division is the part where everything has been turned based on the perspective
        effect (things further away are smaller/vertices are closer to the center of the viewport) */

        this.mapNdcVerticesToScreenCoordinates(verticesData);

        /* Now in 2D screen space. The vertices are finally mapped to the 2D screen coordinates. In this space, we can finally
        see the final position of the vertices. */

        this.displayTrianglesWireframe(trianglesData, verticesData);
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
            new Triangle(this.vertices[0], this.vertices[1], this.vertices[2], Color.RED),
            new Triangle(this.vertices[0], this.vertices[2], this.vertices[3], Color.RED),
            new Triangle(this.vertices[4], this.vertices[5], this.vertices[6], Color.GREEN),
            new Triangle(this.vertices[4], this.vertices[6], this.vertices[7], Color.GREEN),
            new Triangle(this.vertices[0], this.vertices[1], this.vertices[5], Color.BLUE),
            new Triangle(this.vertices[0], this.vertices[5], this.vertices[4], Color.BLUE),
            new Triangle(this.vertices[3], this.vertices[2], this.vertices[6], Color.RED),
            new Triangle(this.vertices[3], this.vertices[6], this.vertices[7], Color.RED),
            new Triangle(this.vertices[0], this.vertices[4], this.vertices[3], Color.GREEN),
            new Triangle(this.vertices[4], this.vertices[3], this.vertices[7], Color.GREEN),
            new Triangle(this.vertices[1], this.vertices[5], this.vertices[6], Color.BLUE),
            new Triangle(this.vertices[1], this.vertices[2], this.vertices[6], Color.BLUE)
        ];
    }
}

export class TriangleMesh extends Mesh {
    constructor(position, size, color) {
        super();
        this.position = position;
        this.size = size;

        this.vertices = [
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y - this.size.y / 2 + 1, this.position.z - this.size.z / 2 - 1)),
        ];

        this.triangles = [
            new Triangle(this.vertices[0], this.vertices[1], this.vertices[2], color),
        ];
    }
}